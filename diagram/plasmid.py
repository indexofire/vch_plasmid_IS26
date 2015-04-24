from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO


genomes = ['0001', '0005', '0007', '0008', '0009', '0010', '0011', '0014', '0015', '0016',
    '0017','0018', '0019', '0020', '0021', '0022', '0024', '0035']
gd_diagram = GenomeDiagram.Diagram("Vibrio cholerae plasmid 2012EL-2176")
max_len = 0

for genome in genomes:
    record = SeqIO.read("../sequences/%s.gb" % genome, "genbank")
    max_len = max(max_len, len(record))
    gd_track_for_features = gd_diagram.new_track(1, greytrack=True, start=0, end=len(record), name=record.name)
    gd_feature_set = gd_track_for_features.new_set()
    for feature in record.features:

        if feature.type != "gene":
            if isinstance(feature.qualifiers.get('product'), list) and 'transposase' in feature.qualifiers.get('product')[0] :

                if len(gd_feature_set) % 2 == 0:
                    color = '#44ACAC'
                else:
                    color = '#EC5503'
                gd_feature_set.add_feature(feature, label=False, sigil="ARROW", color=color,
                    arrowshaft_height=0.25, name=feature.qualifiers.get('product')[0],
                    label_position="middle", label_angle=45)

gd_diagram.draw(format="linear", orientation="landscape", pagesize=(40*cm, 80*cm), fragments=1, start=0, end=max_len)
gd_diagram.write("linear.pdf", "PDF")