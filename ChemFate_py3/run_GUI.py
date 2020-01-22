from tkinter import *
import os
from model_setup import Model_SetUp
from load_data import LoadData
from load_data_nano import load_data
from xlrd import open_workbook

CUR_PATH = os.path.dirname(os.path.abspath(__file__))

def run():
    val1 = get_val1() # chem_type
    chem_type = get_chemType(val1)
    run_option = 1
    bgPercOption2 = 10
    chem_file = get_chemFile(chem_type)
    region_file = CUR_PATH + './Input/Region.xlsx'
    release_file = CUR_PATH + './Input/ChemRelease.xlsx'
    start_date_txt = real_start_date.get()
    end_date_txt = real_end_date.get()
    file_name = outfile_name.get()
    output_file_path = CUR_PATH + './Output/' + file_name + '/'

    if not os.path.exists(output_file_path):
        os.makedirs(output_file_path)

    model = Model_SetUp(start_date_txt, end_date_txt, run_option, bgPercOption2,
                        chem_type, chem_file, region_file, release_file, output_file_path, file_name)
    model.run_model()

    finish_txt = Label(text="ChemFate just finished running, please check the Output folder for results.",
                       font='none 12 bold', bg=bg_col, fg='orange')
    finish_txt.place(x=100, y=530)

    window.after(5000, destroy_widget, finish_txt) # stay for 5 seconds


def verify_input():
    # get the data from input folder
    val1 = v1.get()
    chem_type = get_chemType(val1)
    chem_file = get_chemFile(chem_type)
    region_file = CUR_PATH + './Input/Region.xlsx'
    release_file = CUR_PATH + './Input/ChemRelease.xlsx'
    start_date_temp = '2005 1 1'
    end_date_temp = '2005 12 31'

    if chem_type != 'Nanomaterial':
        load_data_nonNano = LoadData(chem_type, chem_file, region_file, release_file, start_date_temp, end_date_temp)
        chem_params, presence, env, climate, bgConc, release, release_scenario = load_data_nonNano.run_loadData()
    else:
        time, presence, env, climate, bgConc, chem_params, release, release_scenario = load_data(region_file,
                                                                                                release_file,
                                                                                                chem_file,
                                                                                                start_date_temp,
                                                                                                end_date_temp)

    chem_txt4 = Label(text='Chemical:', bg=bg_col).place(x=140, y=170)
    chem_name = StringVar(window, value=chem_params['name'])
    entry1 = Entry(window, textvariable=chem_name, bg='lightgrey', state='disabled')
    entry1.place(x=220, y=170, width=250)

    chem_txt5 = Label(text='Region:', bg=bg_col).place(x=140, y=200)
    region_name = StringVar(window, value=env['name'])
    entry2 = Entry(window, textvariable=region_name, bg='lightgrey', state='disabled')
    entry2.place(x=220, y=200, width=250)

    chem_txt6 = Label(text='Release:', bg=bg_col).place(x=140, y=230)
    release_scenario = StringVar(window, value=release_scenario)
    entry3 = Entry(window, textvariable=release_scenario, bg='lightgrey', state='disabled')
    entry3.place(x=220, y=230, width=250)

    return chem_params, env, release_scenario



def destroy_widget(widget):
    widget.destroy()

def get_val1():
    val = v1.get()
    return val

def get_date(region_file):
    region_workbook = open_workbook(region_file)
    climate = region_workbook.sheet_by_name('Climate')
    start_month = climate.cell_value(rowx=1, colx=0)
    start_day = climate.cell_value(rowx=1, colx=1)
    start_year = climate.cell_value(rowx=1, colx=2)
    end_month = climate.cell_value(rowx=climate.nrows - 1, colx=0)
    end_day = climate.cell_value(rowx=climate.nrows - 1, colx=1)
    end_year = climate.cell_value(rowx=climate.nrows - 1, colx=2)

    start_date = str(int(start_year)) + ' ' + str(int(start_month)) + ' ' + str(int(start_day))
    end_date = str(int(end_year)) + ' ' + str(int(end_month)) + ' ' + str(int(end_day))
    return start_date, end_date

def get_chemFile(chem_type):
    if chem_type == 'NonionizableOrganic':
        chem_file = CUR_PATH + './Input/ChemParam_nonionizableOrganic.xlsx'
    elif chem_type == 'IonizableOrganic':
        chem_file = CUR_PATH + './Input/ChemParam_ionizableOrganic.xlsx'
    elif chem_type == 'Metal':
        chem_file = CUR_PATH + './Input/ChemParam_metal.xlsx'
    else:
        chem_file = CUR_PATH + './Input/ChemParam_nanomaterial.xlsx'
    return chem_file

def get_chemType(val):
    if val == 1:
        chem_type = 'NonionizableOrganic'
    elif val == 2:
        chem_type = 'IonizableOrganic'
    elif val == 3:
        chem_type = 'Metal'
    elif val == 4:
        chem_type = 'Nanomaterial'
    return chem_type

#############################################################################################################
# main
window = Tk()
window.title("ChemFate")
window.geometry("800x620")
bg_col = '#e2f3ff' # '#b9def7'
window.configure(background=bg_col)

# display text
welcome_txt = Label(text="Welcome to ChemFate", font='none 22 bold', fg='royalblue', bg=bg_col)
welcome_txt.place(x=250, y=5) # place in the center

# selection for chemical type
step1_txt = "Step 1. ChemFate contains four fate and transport models: 1) organoFate, for non-ionizable organic chemicals, \n" \
            "         2) ionOFate, for ionizable organic chemicals, 3) metalFate, for metals, 4) nanoFate, for nanomaterials. "
chem_txt1 = Label(text=step1_txt, bg=bg_col).place(x=20, y=50)
step1_txt2 = "            Please select the correct model for your chemical:"
chem_txt2 = Label(text=step1_txt2, bg=bg_col).place(x=20, y=82)

# radio button
v1 = IntVar()
v1.set(1)
radio_button1 = Radiobutton(window, text="organoFate", value=1, variable=v1,
                            command=get_val1, bg=bg_col).place(x=80, y=110)
radio_button2 = Radiobutton(window, text="ionOFate", value=2, variable=v1,
                            command=get_val1, bg=bg_col).place(x=200, y=110)
radio_button3 = Radiobutton(window, text="metalFate", value=3, variable=v1,
                            command=get_val1, bg=bg_col).place(x=320, y=110)
radio_button4 = Radiobutton(window, text="nanoFate", value=4, variable=v1,
                            command=get_val1, bg=bg_col).place(x=440, y=110)

start_date = StringVar()
end_date = StringVar()
# entry for start and end date
step2_txt1 = "Step 2. Please verify the inputs from the files in the 'Input Folder':"
chem_txt3 = Label(text=step2_txt1, bg=bg_col).place(x=20, y=140)

verify_button = Button(window, text="Verify Input", font='none 13 bold', activebackground='cyan', highlightbackground='lightgreen')
verify_button.config(height = 1, width = 12)
verify_button.config(command=verify_input)
verify_button.place(x=500, y=150)


step3_txt1 = "Step 3. Please enter the start date and end date for your model simulation time: \n" \
             "(Note the format: start date: 2005 1 1, end date: 2014 12 31)"

CUR_PATH = os.path.dirname(os.path.abspath(__file__))
region_file = CUR_PATH + './Input/Region.xlsx'
start_date, end_date = get_date(region_file)
chem_txt7 = Label(text=step3_txt1, bg=bg_col).place(x=20, y=270)
chem_txt8 = Label(text="Current date range", bg=bg_col).place(x=190, y=310)
chem_txt9 = Label(text="Start Date:", bg=bg_col).place(x=80, y=340)
chem_txt10 = Label(text="End Date:", bg=bg_col).place(x=80, y=370)
chem_txt11 = Label(text="Desired date range", bg=bg_col).place(x=360, y=310)

start_date = StringVar(window, value=start_date)
entry4 = Entry(window, textvariable=start_date, bg='lightgrey', state='disabled')
entry4.place(x=170, y=340, width=150)

end_date = StringVar(window, value=end_date)
entry4 = Entry(window, textvariable=end_date, bg='lightgrey', state='disabled')
entry4.place(x=170, y=370, width=150)

real_start_date = StringVar()
real_end_date = StringVar()
entry5 = Entry(window, textvariable=real_start_date).place(x=350, y=340, width=150)
entry6 = Entry(window, textvariable=real_end_date).place(x=350, y=370, width=150)

step4_txt1 = "Step 4. Please enter a unique name for your output files (e.g., benzene_SF):"
chem_txt12 = Label(text=step4_txt1, bg=bg_col).place(x=20, y=410)
chem_txt13 = Label(text="Output File Name:", bg=bg_col).place(x=80, y=440)
outfile_name = StringVar()
entry7 = Entry(window, textvariable=outfile_name).place(x=220, y=440, width=300)

step5_txt1 = "Step 5. Ready to run the model:"
chem_txt14 = Label(text=step5_txt1, bg=bg_col).place(x=20, y=480)


# create a button
# background color may not work for Mac, so use highlightbackground
run_button = Button(window, text="Run ChemFate", font='none 14 bold', activebackground='cyan', highlightbackground='lightgreen')
run_button.config(height = 1, width = 15)
run_button.config(command=run)
run_button.place(x=250, y=480)

credit_txt1 = "nanoFate was developed by Dr. Kendra Garner and Dr. Arturo Keller"
credit_txt2 = "organoFate, ionOFate, and metalFate were developed by Dr. Mengya Tao and Dr. Arturo Keller"
credit_txt3 = "Questions: email arturokeller@ucsb.edu; Keller Lab, Bren School, UC Santa Barbara, USA"
chem_txt15 = Label(text=credit_txt1, bg=bg_col, font='none 10 italic').place(x=400, y=555)
chem_txt16 = Label(text=credit_txt2, bg=bg_col, font='none 10 italic').place(x=245, y=573)
chem_txt17 = Label(text=credit_txt3, bg=bg_col, font='none 10 italic').place(x=262, y=592)

# run the main loop
window.mainloop()