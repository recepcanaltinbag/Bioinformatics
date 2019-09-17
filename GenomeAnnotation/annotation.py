#Recep Can Altınbağ
#2019
#Gene annotation visualization of GFF Files with using pyGUI, GFF Parser, Biopython and DnaFeaturesViewer



from dna_features_viewer import GraphicFeature, GraphicRecord
import PySimpleGUI as sg
import os
from BCBio import GFF

cut_number = 15000


in_file = "genome.gff"
in_handle = open(in_file)

list_contigs = []
for rec in GFF.parse(in_handle):
    list_contigs.append(rec.id)
in_handle.close()

list_contigs.reverse()

layout = [ [sg.Text('Pick a Cut Number(between 10000-15000 is ideal)\nIf DNA segment is bigger than cut number, it will divide to parts')],
           [sg.InputText('12000')],
           [sg.Text('Pick a Contig/DNA Segment')], [sg.Listbox(values=list_contigs, size=(50, 20))],
           [sg.Submit(), sg.Text('created by Rca', size = (50,1), justification='right')] ]

window = sg.Window('Annotation Viewer version 0.0').Layout(layout)
button, values = window.Read()

print(values)
cut_number = int(values[0])
print(values[0])
input()
sequence_name = str(values[1][0])
print(sequence_name)

in_file = "genome.gff"
in_handle = open(in_file)
limit_info = dict(
    gff_id = [sequence_name]
)

output_file_name = sequence_name


path = sequence_name
if not os.path.exists(path):
    os.mkdir(path)



for rec in GFF.parse(in_handle, limit_info=limit_info):
    if len(rec.features) > 0:
        features = []
        part = ''
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
        #record.plot(figure_width=5)



        if int(len(rec.seq)) > cut_number:
            element = 0
            while element < int(len(rec.seq)):
                part = str(element)
                zoom_start, zoom_end = element, element+cut_number
                print(element)
                if element + cut_number > len(rec.seq):
                    zoom_end = len(rec.seq) - 1
                    element = len(rec.seq)
                else:
                    element = element + cut_number
                part = 'from' + part + 'to' + str(element)
                output_file_name = path + '/' + sequence_name + part + '.png'
                cropped_record = record.crop((zoom_start,zoom_end))
                ax, _ = cropped_record.plot(figure_width=8)
                ax.figure.savefig(output_file_name, bbox_inches='tight')

        else:
            record.plot(figure_width=5)
            ax, _ = record.plot(figure_width=8)
        #record.plot_translation(ax, (8, 12), fontdict={'weight': 'bold'})

            output_file_name = path + '/' + sequence_name + part + '.png'
            ax.figure.savefig(output_file_name, bbox_inches='tight')

in_handle.close()

popup = "A file is created in working folder named as:  " + output_file_name
sg.Popup(popup)







