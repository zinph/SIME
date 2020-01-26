#-------------------------------------------------------------------------------
# Name:        SIME
# Purpose:     Scaffold-based enumeration method/software to design in-silico macrolide libraries

#
# Author:      zinph
#
# Created:     27/02/2019
# Copyright:   (c) kzphy 2019
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import time
from flask import Flask, render_template, url_for, request
from werkzeug import secure_filename
import os
from SIME import *

app = Flask(__name__)
UPLOAD_FOLDER = 'uploads'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

@app.route('/')
def main_page():
    return render_template('index.html')

@app.route('/collect_data',methods=['POST'])
def collect_data():
    form_data              = request.form
    form_files             = request.files
    macrolide_core         = form_files.get('macrolide_core')
    structural_motifs_file = form_files.get('structural_motifs_file')
    sugars_file            = form_files.get('sugars_file')

    #more options
    max_repeat_motifs      = int(form_data.get("max_repeat_motifs"))
    minimal_sugars         = int(form_data.get("minimal_sugars"))
    library_size           = int(form_data.get("library_size"))
    if form_data.get("enumerate_all_SMs"):
        enumerate_all_SMs  = form_data.get("enumerate_all_SMs")
    else:
        enumerate_all_SMs  = 'no'
    if form_data.get("enumerate_all_sugars"):
        enumerate_all_sugars   = form_data.get("enumerate_all_sugars")
    else:
        enumerate_all_sugars = 'no'
    #enumarate_all          = form_data["enumerate_all"]

    if macrolide_core:
        filename = secure_filename(macrolide_core.filename)
        macrolide_core.save(os.path.join(app.config["UPLOAD_FOLDER"],filename))
        macrolide_core = open(os.path.join(UPLOAD_FOLDER,filename),'r')
    else:
        macrolide_core = None
    if structural_motifs_file:
        filename = secure_filename(structural_motifs_file.filename)
        structural_motifs_file.save(os.path.join(app.config["UPLOAD_FOLDER"],filename))
        structural_motifs_file =open(os.path.join(UPLOAD_FOLDER,filename),'r')
    else:
        structural_motifs_file = None
    if sugars_file:
        filename = secure_filename(sugars_file.filename)
        sugars_file.save(os.path.join(app.config["UPLOAD_FOLDER"],filename))
        sugars_file = open(os.path.join(UPLOAD_FOLDER,filename),'r')
    else:
        sugars_file = None

    sample = SIME(structural_motifs_file, sugars_file, max_repeat_motifs, minimal_sugars, library_size, enumerate_all_SMs, enumerate_all_sugars)
    if macrolide_core == None:
        with open("Data/ery_core.txt", 'r') as f:
            smile = f.readline()
    else:
        smile = macrolide_core.readline()
    start_time = time.time()
    sample.generate_templates_withExtendersNSugars(smile)
    duration = convert_time(time.time()-start_time)
    f'Time Elapsed for Enumeration: {duration}'

    return f'''Time Elapsed for Enumeration: {duration}.
    Your chemical libraries have been generated.
    Please check in LIBRARIES folder. The resulting files for info and smiles should be there.'''
    #return form_data, form_files

def convert_time(second):
    '''
    Convert second to DAY, HOUR, MIN, SEC format. Parameter is second.
    '''
    day = second/86400
    hour = (day - int(day))*24
    minute = (hour - int(hour))*60
    second = round((minute - int(minute))*60,4)
    return(str(int(day)) + ' DAYS: '+ str(int(hour)) + ' HOURS: '+ str(int(minute)) + ' MINUTES: ' + str(second) + ' SECONDS')


if __name__ == '__main__':

    app.jinja_env.auto_reload = True
    app.config['TEMPLATES_AUTO_RELOAD']=True        #forces flask to reload html templates
    app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0     #prevents browsers from caching static files served by flask, such as js code

    app.run(debug=True,use_reloader=True)
    #app.run(host='0.0.0.0',use_reloader=True)
