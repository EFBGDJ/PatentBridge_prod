from flask import Flask
from flask import render_template, url_for,redirect
from flask import request
from flask import session
import pandas as pd
import requests
from Bio.Seq import Seq

def get_token(client_id, client_secret):
    # get credentials from a ENV
    
    access_url = 'https://login.microsoftonline.com/fcb2b37b-5da0-466b-9b83-0014b67a7c78/oauth2/v2.0/token'
    # apiURL = creds['KD_SCORE_DATA_URL']
    # Build the call to Oauth client from creds file
    headers = {'Content-Type': 'application/x-www-form-urlencoded'}
    
    payload = {
                'grant_type' : 'client_credentials',
                'client_id' : client_id,
                'client_secret' : client_secret,
                'scope' : client_id + "/.default"
                }
    response = requests.post(access_url, headers=headers,data=payload)
    # check for a successful response
    if response.status_code == 200:
        accessObj = json.loads(response.content)
        accessToken = accessObj['access_token']
        return (accessToken)
    
    return (None)

def run_vel_blast(sequences, access_token, algo="blastn", collection="nr", evalue=10, nhits=500, word_size=11):
    url="https://velocity.ag/ncbi-blast/services/v1/blast/"
    headers = {
        'accept': 'application/json',
        'Content-Type': 'application/json',
        'authorization': 'Bearer ' + access_token
    }
    user=os.getlogin()
    job=user+"_"+str(datetime.datetime.now()).replace(" ","_")
    param_str="-evalue "+str(float(evalue))+" -num_alignments "+str(nhits)+" -word_size "+str(word_size)
    body={ "userId": user,
           "jobName": job,
           "jobDescription": job+"_"+collection,
           "blastCompletionQueue": "pd-genomics-ncbi-blast-completion-out-srgrod",
           "blastDetails": [ { "algorithm": algo,
                               "parameters": param_str,
                               "collectionName": collection, 
                               "blastFormatters": [ 0,10 ], 
                               "sequences": sequences } ]}
    
    response = requests.post(url, headers=headers, data=json.dumps(body))
    if response.status_code != 200:
        print("Couldn't submit Velocity BLAST job. Error code: "+str(response.status_code))
        return None
    
    Obj = json.loads(response.content)
    batch_id = Obj['blastExecutionBatchId']
    return batch_id


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
