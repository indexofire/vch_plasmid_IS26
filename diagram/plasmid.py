from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO


genomes = ['0001', '0005', '0007', '0008', '0009', '0010', '0011',
    '0014', '0015', '0016', '0017', '0018', '0019', '0020', '0021',
    '0022', '0024', '0035', '0036', '0037', '0039', '0040', '0041',
    '0042', '0043', '0044']

is26 = {
    '0001': ((35928, 36746), (39229, 40048), (151907, 152726)),
    '0005': ((114568, 115387), (118914, 119733), (125526, 126345)),
    '0007': ((156624, 157443), (31272, 32091), (44315, 45134), (136266, 137085),
        (43071,43890), (38573, 39392), (145083, 145902), (33104, 33923), (149429, 150248)),
    '0008': ((5389, 6208), (4347, 5166), (1, 820)),
    '0009': ((362, 1181), (12947, 13766),(11703, 12522), (6540, 7359), (2194, 3013)),
    '0010': ((20337, 21156), (141275, 142094), (19093, 19912)),
    '0011': ((12304, 13123), (19277, 20096), (1460, 2279)),
    '0014': ((117640, 118459), (110574, 111393), (112434, 113253)),
    '0015': ((4517, 5336), (19042, 19861), (150016, 150835), (25598, 26417), (23738, 24557),
        (179365, 180184), (165323, 166142), (165323, 166142)),
    '0016': ((157299, 158118), (20336, 21155), (19093, 19912), (127230, 128049), (113483,114302)),
    '0017': ((22025, 22844), (131331, 132150), (20781, 21600)),
    '0018': ((25181, 26000), (23937, 24756), (18568, 19387)),
    '0019': ((20337, 21156), (126141, 126959), (19093, 19912)),
    '0020': ((20337, 21156), (125881, 126700), (19094, 19913)),
    '0021': ((20337, 21156), (123015, 123834), (19093, 19912)),
    '0022': ((25979, 26798), (24731, 25550)),
    '0024': ((17257, 18076),),
    '0035': ((19220, 20039), (137161, 137980), (13976, 18795)),
    '0036': ((176068, 176887), (160010, 160829), (163663, 164482), (167558, 168377), (60970, 61789), (211083, 211906), (219904,220723), (32816, 33635)),
    '0037': ((22056, 22875), (87666, 88485), (20812, 21631), (100736, 101555)),
    '0039': ((1, 820), (53526, 54345), (147857, 148643)),
    '0040': ((143281, 144100), (45916, 46734)),
    '0041': ((223747, 224566), (207689, 208508), (211342, 212161), (215237, 216056), (60970, 61789), (258766, 259585), (267583, 268402), (32816, 33635)),
    '0042': ((176079, 176898), (160021, 160840), (163674, 1644493), (211098, 211917), (219915, 220734), (32816, 33635)),
    '0043': ((49405, 50224), (36081, 36900), (40427, 41246)),
    '0044': ((20337, 21156), (106856, 107675)),
}

is26_truncated = {
    '0007': ((30894, 31271),),
}

gd_diagram = GenomeDiagram.Diagram("Vibrio cholerae plasmid 2012EL-2176")
max_len = 0

for genome in genomes:
    record = SeqIO.read("../sequences/%s.gb" % genome, "genbank")
    max_len = max(max_len, len(record))
    gd_track_for_features = gd_diagram.new_track(1, greytrack=True, start=0,
        end=len(record), name=record.name)
    gd_feature_set = gd_track_for_features.new_set()
    for feature in record.features:

        if feature.type != "gene":
            if isinstance(feature.qualifiers.get('product'), list) and 'transposase' in feature.qualifiers.get('product')[0]:

                #if len(gd_feature_set) % 2 == 0:
                #    color = '#44ACAC'
                #else:
                #    color = '#EC5503'
                color = '#33ACAC'
                name = feature.qualifiers.get('product')[0]

                for is_pos in is26[genome]:
                    if feature.location.start > is_pos[0] and feature.location.end < is_pos[1]:
                        color = '#DBCC36'
                        name = 'IS26'

                #if is26_truncated.setdefault(genome) is not None:
                #    for is_truncated_pos in is26_truncated[genome]:
                #        print feature.location.start, feature.location.end
                #        if feature.location.start == is_truncated_pos[0] and feature.location.end ==is_truncated_pos[1]:
                #            print is_truncated_pos
                #            color = '#66EC33'
                #            name = 'is26_truncated'
                #print feature.location.start, feature.location.end
                gd_feature_set.add_feature(feature, label=False, sigil="ARROW", color=color,
                    arrowshaft_height=0.4, arrowhead_length=0.3, name=name,
                    label_position="middle", label_angle=45)

            elif 'transposase' in str(feature.qualifiers.get('note')):
                color = '#693ABB'
                name = feature.qualifiers.get('product')

                for is_pos in is26[genome]:
                    if feature.location.start > is_pos[0] and feature.location.end < is_pos[1]:
                        color = '#DBCC36'
                        name = 'IS26'

                gd_feature_set.add_feature(feature, label=False, sigil="ARROW", color=color,
                    arrowshaft_height=0.4, arrowhead_length=0.3, name=name,
                    label_position="middle", label_angle=45)

gd_diagram.draw(format="linear", orientation="landscape", pagesize=(40*cm, 160*cm), fragments=1, start=0, end=max_len)
gd_diagram.write("linear.pdf", "PDF")
