import os
from collections import defaultdict

# read GFF file and structure data
contigs = {}

with open('GCF_003789085.1_ASM378908v1_genomic.gff', 'r') as f_input:
    for line in f_input:
        if line.startswith('#') or not line.strip():
            continue

        parts = line.strip().split('\t')
        contig_id = parts[0]
        feature_type = parts[2]
        start = int(parts[3])
        end = int(parts[4])
        strand = parts[6]
        metadata = parts[8]
        metadata_split = metadata.split(';')
        metadata_dict = {item.split('=')[0]: item.split('=')[1] for item in metadata_split if '=' in item}

        feature = {
            'type': feature_type,
            'start': start,
            'end': end,
            'strand': strand,
            'metadata': metadata_dict,
        }

        if contig_id not in contigs:
            contigs[contig_id] = []
        contigs[contig_id].append(feature)

# load the IDs of interest from a list
input_filename = 'ids_lnctar_output.txt'
ids_interes = set()
with open(input_filename, 'r') as f_ids:
    for line in f_ids:
        ids_interes.add(line.strip())

print('Número de lncRNAs cargados:', len(ids_interes))

# load the list of DE mRNAs derived from de DEGs analysis
de_mrnas = set()
with open('de_mRNAs.txt', 'r') as f_de:
    for line in f_de:
        de_mrnas.add(line.strip())

print('Número de mRNAs DE:', len(de_mrnas))

# find the contigs where the IDs of interest are present
contigs_interes = defaultdict(list)

for contig_id, features in contigs.items():
    for feature in features:
        if feature['type'] == 'lnc_RNA':
            metadata_dict = feature['metadata']
            gene_id = metadata_dict.get('ID', 'Unknown')
            if gene_id.startswith('rna-'):
                gene_id = gene_id[4:]  # remove the 'rna-' prefix
            if gene_id in ids_interes:
                contigs_interes[gene_id].append(contig_id)

# this is the main function (search for multiple neighboring mRNAs, 5 in this case)
def find_multiple_mRNA_neighbors(features, lncRNA_index, num_neighbors=5):
    """
    Search multiple neighboring mRNAs for a given lncRNA.
    Only consider the first mRNA after a gene entry as a new gene,
    ignoring mRNA variants.
    """
    upstream_mRNAs = []
    downstream_mRNAs = []
    lncRNA = features[lncRNA_index]

    last_gene_seen = None
    for i, feature in enumerate(features):
        if feature['type'] == 'gene':
            last_gene_seen = i  # mark the last gene position
        elif feature['type'] == 'mRNA' and last_gene_seen is not None:
            if feature['start'] > lncRNA['end']:  # downstream mRNA
                downstream_mRNAs.append(i)
                last_gene_seen = None  # reset after considering the mRNA
            elif feature['end'] < lncRNA['start']:  # upstream mRNA
                upstream_mRNAs.append(i)
                last_gene_seen = None  # reset after considering the mRNA

    # limit to the desired number of neighbors (5 for this time)
    upstream_mRNAs = sorted(upstream_mRNAs, key=lambda x: features[x]['end'], reverse=True)[:num_neighbors]
    downstream_mRNAs = sorted(downstream_mRNAs, key=lambda x: features[x]['start'])[:num_neighbors]

    return upstream_mRNAs, downstream_mRNAs

# calculate the distance between lncRNA and mRNA in Kb
def calculate_distance(lncRNA, mRNA):
    if mRNA['end'] < lncRNA['start']:
        return (lncRNA['start'] - mRNA['end']) / 1000
    elif mRNA['start'] > lncRNA['end']:
        return (mRNA['start'] - lncRNA['end']) / 1000
    else:
        return 0  # overlapping
    

# output filename
output_filename = os.path.splitext(input_filename)[0] + "_rnacontigprofiler_output.txt"


# open the output file for writing
with open(output_filename, 'w') as f_output:
    # write the header
    f_output.write("lncRNA_ID\tstrand\tcontig\tmRNA_ID\tmRNA_Position\tDE_status\tdirection\tmRNA_strand\tproduct\tdistance_Kb\n")

    # search for neighboring mRNAs for the IDs of interest and check if they are DE
    for gene_id in ids_interes:
        for contig_id in contigs_interes[gene_id]:
            features = contigs[contig_id]
            for i, feature in enumerate(features):
                if feature['type'] == 'lnc_RNA' and feature['metadata'].get('ID', '').endswith(gene_id):
                    upstream_mRNAs, downstream_mRNAs = find_multiple_mRNA_neighbors(features, i)

                    # process upstream mRNAs
                    if upstream_mRNAs:
                        for idx, upstream_mRNA_idx in enumerate(upstream_mRNAs):
                            upstream_mRNA = features[upstream_mRNA_idx]
                            upstream_mRNA_id = upstream_mRNA['metadata'].get('ID', 'Unknown')
                            upstream_mRNA_product = upstream_mRNA['metadata'].get('product', 'Unknown')
                            if upstream_mRNA_id.startswith('rna-'):
                                upstream_mRNA_id = upstream_mRNA_id[4:]  # Remove the 'rna-' prefix
                            distance_upstream = calculate_distance(feature, upstream_mRNA)
                            de_status_upstream = "DE" if upstream_mRNA_id in de_mrnas else "not DE"
                            f_output.write(f"{gene_id}\t{feature['strand']}\t{contig_id}\t{upstream_mRNA_id}\t{idx + 1}_U\t{de_status_upstream}\tupstream\t{upstream_mRNA['strand']}\t{upstream_mRNA_product}\t{distance_upstream:.3f}\n")

                    # process downstream mRNAs
                    if downstream_mRNAs:
                        for idx, downstream_mRNA_idx in enumerate(downstream_mRNAs):
                            downstream_mRNA = features[downstream_mRNA_idx]
                            downstream_mRNA_id = downstream_mRNA['metadata'].get('ID', 'Unknown')
                            downstream_mRNA_product = downstream_mRNA['metadata'].get('product', 'Unknown')
                            if downstream_mRNA_id.startswith('rna-'):
                                downstream_mRNA_id = downstream_mRNA_id[4:]  # Remove the 'rna-' prefix
                            distance_downstream = calculate_distance(feature, downstream_mRNA)
                            de_status_downstream = "DE" if downstream_mRNA_id in de_mrnas else "not DE"
                            f_output.write(f"{gene_id}\t{feature['strand']}\t{contig_id}\t{downstream_mRNA_id}\t{idx + 1}_D\t{de_status_downstream}\tdownstream\t{downstream_mRNA['strand']}\t{downstream_mRNA_product}\t{distance_downstream:.3f}\n")

print(f"Results have been written to {output_filename}")
