#Recep Can Altınbağ
#2019
#Gene annotation visualization of GFF Files with using pyGUI, GFF Parser, Biopython and DnaFeaturesViewer

from dna_features_viewer import GraphicFeature, GraphicRecord
import PySimpleGUI as sg
import pprint
from BCBio.GFF import GFFExaminer
from BCBio import GFF


#Open File
in_file = "genome.gff"
in_handle = open(in_file)

list_contigs = []
for rec in GFF.parse(in_handle):
    list_contigs.append(rec.id)
in_handle.close()

list_contigs.reverse()

#Very Basic User Interface with pyGUI
layout = [ [sg.Text('Pick a Contig/DNA Segment')], [sg.Listbox(values=list_contigs, size=(50, 20))],
           [sg.Submit(), sg.Text('created by Rca', size = (50,1), justification='right')] ]

window = sg.Window('Annotation Viewer version 0.0').Layout(layout)
button, values = window.Read()


print(values[0][0])
sequence_name = str(values[0][0])
print(sequence_name)

in_file = "genome.gff"
in_handle = open(in_file)
limit_info = dict(
    gff_id = [sequence_name]
)


#Parsing and visualization, you can select the colors you want  
output_file_name = sequence_name + '.png'
for rec in GFF.parse(in_handle, limit_info=limit_info):
    if len(rec.features) > 0:
        features = []
        for feature in rec.features:

            set_color = '#ccccff'
            print(feature.type)
            if str(feature.qualifiers["product"][0]) == 'hypothetical protein':
                set_color = '#ffcccc'
            if feature.type == 'tRNA':
                set_color = "#ff0000"

            print(str(feature.qualifiers["product"]),1)
            features.append(GraphicFeature(start=feature.location.start,
                           end=feature.location.end,
                           strand=feature.location.strand,
                           color=set_color, label=str(feature.qualifiers["product"][0])))

        record = GraphicRecord(sequence=str(rec.seq), features=features)
        record.plot(figure_width=5)

        ax, _ = record.plot(figure_width=8)
        #record.plot_translation(ax, (8, 12), fontdict={'weight': 'bold'})

        ax.figure.savefig(output_file_name, bbox_inches='tight')

in_handle.close()

popup = "A file is created in working folder named as:  " + output_file_name
sg.Popup(popup)



