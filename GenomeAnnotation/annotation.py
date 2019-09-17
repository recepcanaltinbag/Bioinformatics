#Recep Can Altınbağ
#2019
#Gene annotation visualization of GFF Files with using pyGUI, GFF Parser, Biopython and DnaFeaturesViewer
#pip install dna_features_viewer
#pip install bokeh pandas
#pip install biopython
#pip install PySimpleGUI
#pip install bcbio-gff

from dna_features_viewer import GraphicFeature, GraphicRecord
import PySimpleGUI as sg
import os
from BCBio import GFF

cut_number = 15000

start_the_program = False

while start_the_program == False:
    layoutfile = [[sg.Text('Write the directory of your .gff file example: Desktop/genome.gff')],
                  [sg.InputText('genome.gff')],
                  [sg.Submit()]]
    windowfile = sg.Window('Annotation Viewer version 0.0 File').Layout(layoutfile)
    button, values = windowfile.Read()
    in_file = values[0]
    print(button)
    if button == 'Submit':
        try:
            in_handle = open(in_file)
            start_the_program = True
            windowfile.close()
        except:
            popup = "Enter valid directory"
            sg.Popup(popup)
    else:
        exit(1)



list_contigs = []
for rec in GFF.parse(in_handle):
    list_contigs.append(rec.id)
in_handle.close()

list_contigs.reverse()
layout = [
    [sg.Text('Pick a Cut Number(between 10000-15000 is ideal)'
             '\nIf DNA segment is bigger than cut number, it will divide to parts'
             '{\nor you can write a specific portion of the contig/sequence')],
               [sg.InputText('12000'),sg.Checkbox('Specific part of the DNA (ex: 12345-12678)')],
               [sg.Text('Pick a Contig/DNA Segment')], [sg.Listbox(values=list_contigs, size=(50, 20))],
               [sg.Submit(), sg.Text('created by Rca', size = (50,1), justification='right')] ]

window = sg.Window('Annotation Viewer version 0.0').Layout(layout)


terminate = False

while terminate == False:

    button, values = window.Read()
    if values[0]==None or values[2]==None:
        exit(1)

    print(values)

    portion_dna = values[1]
    if portion_dna == False:
        cut_number = int(values[0])
    else:
        start_char = int(values[0].split('-')[0])
        end_char = int(values[0].split('-')[1])
        print(start_char)
        print(end_char)
    sequence_name = str(values[2][0])
    print(sequence_name)


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


            print(float(len(rec.seq)/cut_number))

            if portion_dna == True:
                print('yey')
                if end_char<=len(rec.seq):
                    print('yey')
                    zoom_start, zoom_end = start_char, end_char
                    part = str(zoom_start) + 'to' + str(zoom_end)
                    output_file_name = path + '/' + sequence_name + part + '.png'
                    cropped_record = record.crop((zoom_start, zoom_end))
                    ax, _ = cropped_record.plot(figure_width=8)
                    ax.figure.savefig(output_file_name, bbox_inches='tight')

                else:
                    popup = "Enter a portion less than segment length"
                    sg.Popup(popup)
                    continue


            elif int(len(rec.seq)) > cut_number and float(len(rec.seq)/cut_number) > 1.3:
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
                        if float(len(rec.seq)/element) < 1.3:
                            zoom_end = len(rec.seq) - 1
                            element = len(rec.seq)

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










