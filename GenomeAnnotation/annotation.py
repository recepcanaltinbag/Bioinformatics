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

cut_number = 12000
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
             '{\nor you can write a specific portion of the contig/sequence (ex: 12345-12678)')],
               [sg.InputText('12000')],
               [sg.Text('Pick a Contig/DNA Segment')], [sg.Listbox(values=list_contigs, size=(50, 20))],
               [sg.Submit(), sg.Text('created by Rca', size = (50,1), justification='right')] ]

window = sg.Window('Annotation Viewer version 0.0').Layout(layout)


terminate = False

def smart_cutter_function(start,array):

    min = 100000
    value = 0
    final_start = 0
    for element in array:
        if abs(start-element) < min:
            min = abs(start-element)
            final_start = element

    return final_start


while terminate == False:

    smart_cutter = []
    dna_portion = False
    button, values = window.Read()

    if values[0]==None or values[1]==None:
        exit(1)

    if len(values[0].split('-')) <= 1:
        cut_number = int(values[0])
    else:
        start_char = int(values[0].split('-')[0])
        end_char = int(values[0].split('-')[1])
        dna_portion = True
        print(start_char)
        print(end_char)
    sequence_name = str(values[1][0])
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
                my_label = ''
                set_color = '#ccccff'
                print(feature.type)
                if "product" in feature.qualifiers:
                    if str(feature.qualifiers["product"][0]) == 'hypothetical protein':
                        set_color = '#ffcccc'
                    my_label = str(feature.qualifiers["product"][0])
                elif "description" in feature.qualifiers:
                    set_color = '#ffcccc'
                    my_label = str(feature.qualifiers["description"][0])
                elif "ID" in feature.qualifiers:
                    set_color = '#ffcccc'
                    my_label = str(feature.qualifiers["ID"][0])
                else:
                    print(feature.qualifiers)
                    my_label = str(feature.type)

                if feature.type == 'tRNA':
                    set_color = "#ff0000"

                #print(str(feature.qualifiers["product"]),1)

                smart_cutter.append(feature.location.start)
                smart_cutter.append(feature.location.end)

                print(feature.location.start)
                print(feature.location.end)
                features.append(GraphicFeature(start=feature.location.start,
                               end=feature.location.end,
                               strand=feature.location.strand,
                               color=set_color, label=my_label))

            record = GraphicRecord(sequence=str(rec.seq), features=features)
            #record.plot(figure_width=5)


            print(float(len(rec.seq)/cut_number))


            element = 0
            if int(len(rec.seq)) > cut_number and float((len(rec.seq)-element)/cut_number) > 1.3:

                if dna_portion == True:
                    element = smart_cutter_function(start_char,smart_cutter)
                    lentgh_of = smart_cutter_function(end_char,smart_cutter)
                else:
                    element = 0
                    lentgh_of = int(len(rec.seq))

                while element < lentgh_of:
                    part = str(element)
                    zoom_start, zoom_end = element, smart_cutter_function(element + cut_number, smart_cutter)
                    print(element,lentgh_of,'error')
                    if element + cut_number > lentgh_of:
                        zoom_end = lentgh_of - 1
                        element = lentgh_of
                    else:
                        element = smart_cutter_function(element + cut_number, smart_cutter)
                        if float((len(rec.seq)-element)/cut_number) < 1.3:
                            zoom_end = len(rec.seq) - 1
                            element = len(rec.seq)

                    print(zoom_start,zoom_end,"borders")
                    part = 'from' + part + 'to' + str(element)
                    output_file_name = path + '/' + sequence_name + part + '.png'
                    cropped_record = record.crop((zoom_start,zoom_end))
                    ax, _ = cropped_record.plot(figure_width=8)
                    ax.figure.savefig(output_file_name, bbox_inches='tight')

            else:
                print('normal')
                record.plot(figure_width=5)
                ax, _ = record.plot(figure_width=8)
            #record.plot_translation(ax, (8, 12), fontdict={'weight': 'bold'})

                output_file_name = path + '/' + sequence_name + part + '.png'
                ax.figure.savefig(output_file_name, bbox_inches='tight')

    in_handle.close()
    popup = "A file is created in working folder named as:  " + output_file_name
    sg.Popup(popup)
















