from flask import Flask
from flask import render_template, url_for,redirect
from flask import request
from flask import session
import pandas as pd
import requests
import json
import os
import datetime
import time
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
    user= session$user
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

def get_exec_results(batch_id, access_token, poll_time):
    url="https://velocity.ag/ncbi-blast/services/v1/blast/blast-execution-batch/"
    headers = {'accept': 'application/json', 'authorization': 'Bearer ' + access_token}
    #url_details=requests.get(url+str(batch_id),headers=headers)
    err_count = 0
    MAX_ERRS = 3
    while True:
        time.sleep(poll_time)
        response = requests.get(url+str(batch_id),headers=headers)
        # print(response.text)
        if response.ok:
            response_json = response.json()
            response_set = set(map(lambda x: x["status"], response_json['blastExecutionDetails']))
            print(response_set)
            if response_set <= {'Archive', 'Error'}:
                return response_json
        else:
            print(response.text)
            err_count += 1
            if err_count > MAX_ERRS:
                raise ValueError("Error count greater than max errs, stopping")

def get_result(batch_id, access_token, poll_time):
    url="https://velocity.ag/ncbi-blast/services/v1/blast/result?resultUrl="
    headers = {'accept': 'application/json', 'authorization': 'Bearer ' + access_token}
    exec_res=get_exec_results(batch_id, access_token, poll_time)
    
    #results = []
    warnings=[]
    for result_detail in exec_res['blastExecutionDetails']:
        if result_detail["status"] == 'Archive':
            result_url = result_detail['blastFormatResults'][1]['resultUrl']
            
            response = requests.get(url+result_url,headers=headers)
            return response
            #blast_out = list(filter(lambda x: len(x.rstrip()) > 0, response.text.split("\n")))
            #results.append(blast_out)
        else:
            #warnings.warn("Sequence %s has non-archive status" % result_detail["sequenceName"])
            warnings.append("Sequence %s has non-archive status:  %s"
                                        % (result_detail["sequenceName"], result_detail["errorText"]))
            
    return warnings

def get_blast_res(batch_id, access_token):
    try:
        
        headers = {'accept': 'application/json', 'authorization': 'Bearer ' + access_token}
        url='https://velocity.ag/ncbi-blast/services/v1/blast/blast-execution-batch/'+str(batch_id)+'/result/csv'
        res=requests.get(url, headers=headers)
        
        return res
    except:
        print("Couldn't find blast results")


client_id = '6212df18-f0d8-49e6-a8fc-ea98aae348ad'
client_secret = 'FsA8Q~DnvHjxtmEKF4K_5vYaJLVpx6j4.0BpTc2K'

# Generating token by calling the function get_token # 
token=get_token(client_id, client_secret)
print(token)
	
	
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
    n = session.get('user_name',None)
    if(n[0]=='W'):
        seq_list= wipo_file[wipo_file['patent_number']== n]['Protein'].to_list()
    else:
        seq_list=uspto_file[uspto_file['patent_number']==n]['Protein'].to_list()
    seq_list.append("choose all")
 
    print(seq_list)
    return render_template('item.html',list_new = seq_list, user = n)


@app.route('/display')
def display():
    it = session.get('chosen',None)
    sequence=[{"name":"test",
               "text":it}]
    print(sequence)
    bid=run_vel_blast(sequence,token,algo="blastp",collection="IC_TOXIN-All",
                      evalue=10, nhits=250, word_size=3)
    print(bid)
    res=get_result(bid,token,30)
    print(res)
    with open("blas_res5.csv","w") as file:
        file.write(res.text)
    res_file=pd.read_csv("blas_res5.csv",header=None) 
    res_file.to_csv('blas_res5.csv', header=['qacc', 'sacc', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'], index=False)
    df = pd.read_csv("blas_res5.csv")      
    os.remove("blas_res5.csv") 
    return render_template('results.html', seq = it, column_names = df.columns.values,row_data = list(df.values.tolist()),zip = zip)

	
	

if __name__ == '__main__':
	app.run()
