from flask import Flask
from flask import render_template, url_for,redirect
from flask import request
from flask import session
import pandas as pd
import requests
from Bio.Seq import Seq

app = Flask(__name__)
app.secret_key = "hello"

wipo_file = pd.read_csv("fetched_wipo_all_final.csv")
uspto_file = pd.read_csv("fetched_uspto_all_final.csv")
wipo_list = wipo_file['patent_number'].to_list()
uspto_list = uspto_file['patent_number'].to_list()

@app.route('/',methods = ['POST','GET'])
def index():
    list = wipo_list + uspto_list
    new_list = []
    for items in list:
        if items in new_list:
            continue
        else:
            new_list.append(items)
	#list = ['Drinks','Vegetables','Fruits']
    if request.method == 'POST':
        n = request.form['n1']
        session['user_name'] = n
        return redirect(url_for('name'))	
    return render_template('list.html', list = new_list)

@app.route('/name',methods = ['POST','GET'])
def name():
    if request.method == 'POST':
        it = request.form['n2']
        session['chosen'] = it
        return  redirect(url_for('display'))
    seq_list = ['Choose all']
    n = session.get('user_name',None)
    if(n[0]=='W'):
        seq_list= wipo_file[wipo_file['patent_number']== n]['Protein'].to_list()
    else:
        seq_list=uspto_file[uspto_file['patent_number']==n]['Protein'].to_list()
    print(seq_list)
    return render_template('item.html',list_new = seq_list, user = n)


@app.route('/display')
def display():
	return 'Great'

if __name__ == '__main__':
	app.run()
