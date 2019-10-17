from Bio import SeqIO
from dna_features_viewer import GraphicFeature, GraphicRecord
import PySimpleGUI as sg
import os


cut_number = 12000
start_the_program = False
list_contigs = []
list_feature_types = []
in_file = ''


#
# layoutfile =  [[sg.ColorChooserButton('Color')],[sg.Submit()]]
# windowfile = sg.Window('Annotation Viewer version 0.1 File').Layout(layoutfile)
# button, values = windowfile.Read()
# in_file = values['Color']
# print(in_file)



def color_picking(list_feature_types):
    my_color_dict = {}
    i = 0
    color_array = ['#BDB76B','#00CED1','#DEB887','#87CEEB','#F4A460','#FF8C00','#DA70D6','#008080','#000080','#8B0000','#FFB6C1','#EE82EE','#8B008B']

    my_array = sorted(list_feature_types)

    for feature in my_array:
        my_color_dict[feature] = color_array[i]
        i = i + 1
        if i > len(list_feature_types)-1:
            i = 0

    return my_color_dict


my_color_dict = {}

while start_the_program == False:
    layoutfile = [[sg.Text('Write the directory of your .gff file example: Desktop/genome.gbk')],
                  [sg.InputText('annot.gbk')],
                  [sg.Submit()]]
    windowfile = sg.Window('Annotation Viewer version 0.1 File').Layout(layoutfile)
    button, values = windowfile.Read()
    in_file = values[0]
    print(button)
    if button == 'Submit':
        try:
            with open(in_file, "rU") as input_handle:
                for record in SeqIO.parse(input_handle, "genbank"):
                    list_contigs.append(record.id)
                    for feature in record.features:
                        list_feature_types.append(feature.type)
            start_the_program = True

            windowfile.close()
        except:
            popup = "Enter valid directory"
            sg.Popup(popup)
    else:
        exit(1)

list_feature_types = set(list_feature_types)
my_color_dict = color_picking(list_feature_types)
print(my_color_dict)
print(list_feature_types)
print("yey")


list_contigs.reverse()
layout = [
    [sg.Text('Pick a Cut Number(between 10000-15000 is ideal)'
             '\nIf DNA segment is bigger than cut number, it will divide to parts'
             '{\nor you can write a specific portion of the contig/sequence (ex: 12345-12678)')],
               [sg.InputText('12000')],
               [sg.Text('Pick a Contig/DNA Segment')], [sg.Listbox(values=list_contigs, size=(50, 20))],
               [sg.Submit(), sg.Text('created by Rca', size = (50,1), justification='right')] ]

window = sg.Window('Annotation Viewer version 0.1').Layout(layout)
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
    ID = str(values[1][0])
    print(ID)


    output_file_name = ID
    path = ID
    if not os.path.exists(path):
        os.mkdir(path)

    with open(in_file, "rU") as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            if record.id == ID:
                part = ''
                features = []
                feature = 0
                while feature < len(record.features):
                    my_label = ''
                    set_color = '#ccccff'
                    increase_feature = False

                    print('feature',feature)

                    if "source" == record.features[feature].type:
                        feature = feature + 1
                        continue

                    if "gene" == record.features[feature].type:
                        gene_location = record.features[feature].location
                        cds_location = record.features[feature + 1].location
                        strand = record.features[feature].location.strand

                        if "product" in record.features[feature + 1].qualifiers and cds_location == gene_location:
                            start_point = cds_location.start
                            end_point = cds_location.end
                            strand = strand
                            my_label = record.features[feature + 1].qualifiers["product"][0]
                            increase_feature = True


                    else:
                        gene_location = record.features[feature].location
                        strand = record.features[feature].location.strand
                        start_point = gene_location.start
                        end_point = gene_location.end
                        my_label = record.features[feature].type

                    smart_cutter.append(start_point)
                    smart_cutter.append(end_point)
                    set_color = my_color_dict[record.features[feature].type]
                    if my_label == 'hypothetical protein':
                        set_color = '#90EE90'

                    if increase_feature:
                        feature = feature + 2
                    else:
                        feature = feature + 1


                    features.append(GraphicFeature(start=start_point,
                                                       end=end_point,
                                                       strand=strand,
                                                       color=set_color, label=my_label))

                graph_record = GraphicRecord(sequence=str(record.seq), features=features)
                # record.plot(figure_width=5)

                print(float(len(record.seq) / cut_number))


                element = 0
                if int(len(record.seq)) > cut_number and float((len(record.seq)-element)/cut_number) > 1.3:

                    if dna_portion == True:
                        element = smart_cutter_function(start_char,smart_cutter)
                        lentgh_of = smart_cutter_function(end_char,smart_cutter)
                    else:
                        element = 0
                        lentgh_of = int(len(record.seq))

                    while element < lentgh_of:
                        part = str(element)
                        zoom_start, zoom_end = element, smart_cutter_function(element + cut_number, smart_cutter)
                        print(element,lentgh_of,'error')
                        if element + cut_number > lentgh_of:
                            zoom_end = lentgh_of - 1
                            element = lentgh_of
                        else:
                            element = smart_cutter_function(element + cut_number, smart_cutter)
                            if float((len(record.seq)-element)/cut_number) < 1.3:
                                zoom_end = len(record.seq) - 1
                                element = len(record.seq)

                        print(zoom_start,zoom_end,"borders")
                        part = 'from' + part + 'to' + str(element)
                        output_file_name = path + '/' + ID + part + '.png'
                        cropped_record = graph_record.crop((zoom_start,zoom_end))

                        ax, _ = cropped_record.plot(figure_width=8)
                        ax.figure.savefig(output_file_name, bbox_inches='tight')

                else:
                    print('normal')
                    graph_record.plot(figure_width=5)
                    ax, _ = graph_record.plot(figure_width=8)
                #record.plot_translation(ax, (8, 12), fontdict={'weight': 'bold'})

                    output_file_name = path + '/' + ID + part + '.png'
                    ax.figure.savefig(output_file_name, bbox_inches='tight')


    popup = "A file is created in working folder named as:  " + output_file_name
    sg.Popup(popup)







