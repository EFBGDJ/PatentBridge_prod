from flask import Flask,g
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
from flask import make_response

app = Flask(__name__)
app.secret_key = "hello"

def alignment(seq1,seq2):
    file = open("alignment.txt","w")
    aligner = Align.PairwiseAligner(match_score = 1.0)
    score = aligner.score(seq1,seq2)
    alignment = aligner.align(seq1,seq2)
    file.write(str(aligner))
    file.write(f"The selected sequence: {seq1}\n")
    file.write(f"The manually pasted sequence: {seq2}\n")
    file.write(f"The alignment score is {score}\n")
    for var in alignment:
        var = str(var)
        file.write(var)
 
 
# -------- Function to process the entered nucleotide sequence----------- #
def process_string(str):
    val = ''.join([i for i in str if not i.isdigit()])
    val = ''.join(val.splitlines())
    val = val.replace(" ","")
    return val
 
 
# ----------Function to process the entered protein sequence------------#
def process_protein(protein_string):
    dicts = {'Ala': 'a', 'Asx':'b','Cys':'c','Asp':'d','Glu':'e','Phe':'f','Gly':'g','His':'h','Ile':'i','Lys':'k','Leu':'l','Met':'m','Asn':'n','Pro':'p','Gln':'q','Arg':'r','Ser':'s','Thr':'t','Sec':'u','Val':'v','Trp':'w','Xaa':'x','Tyr':'y','Glx':'z'}
    protein_listi = ''.join([i for i in protein_string if not i.isdigit()])
    protein_listi = ''.join(protein_listi.splitlines())
    protein_listi = protein_listi.split(" ")
    c = protein_listi.count('')
    for i in range(c):
        protein_listi.remove('')
    protein_string = ""
    for amino in protein_listi:
        if amino in dicts.keys():
            protein_string += dicts[amino]
    protein_string = protein_string.upper()
    return protein_string
    
    return protein_listi
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
    user= "EFBDJ"
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
        return redirect(url_for('original_sequence'))	
    return render_template('list.html', list = new_list)
 
 
#--- App decorator to fetch the original sequences---#
@app.route('/original_sequence', methods = ['POST','GET'])
def original_sequence():
    global sequences
    if request.method == 'POST':
        if request.form.get('action1') == 'Download sequences':
            return redirect(url_for('download_original'))
        if request.form.get('action2') == 'Translate':
            return redirect(url_for('translate'))
 
    n = session.get('user_name',None)
    if(n[0]=='W'):
        #if n in wipo_file['patent_number']:    
            #sequences = wipo_file[wipo_file['patent_number']== n]['Protein'].to_list()
        #else:
            #result = get.get_seq_wipo(n)
            #print(result)
            #sequences = result['Protein'].to_list()
            #seq_list = result['Protein'].to_list()
        sequences_df = wipo_file[['patent_number','Sequence']]
        sequences_df = sequences_df[sequences_df['patent_number'] == n] 
        sequences = sequences_df['Sequence'].to_list()
    else:
        sequences_df = uspto_file[['patent_number','Sequence']]
        sequences_df = sequences_df[sequences_df['patent_number'] == n] 
        sequences = sequences_df['Sequence'].to_list()
    return render_template('sequence.html', sequence = sequences)
 
 
@app.route('/download_original', methods = ['POST','GET'])
def download_original():
    print(sequences)
    n = session.get('user_name',None)
    print(n)
    index = list(range(1,len(sequences)+1)) 
    fasta_file_string = ",".join(sequences)
    print(fasta_file_string)
        # Replacing the string with new line and > character for FASTA format #
    fasta_file = fasta_file_string.replace(",","\n>\n")
    print(fasta_file)
        # Adding > for the first sequence
    fasta_file = ">\n" + fasta_file
    i = 1
    fasta = ""
    for line in fasta_file:
        if(line == ">"):
            fasta += "\n"
            id = n + "_" + str(i)
        # To give numbers to each sequence #
            new_line = line + id
            fasta += new_line
            i += 1
        else:
            fasta += line
    response = make_response(fasta)
    cd = 'attachment; filename= sequence.fasta'
    response.headers['Content-Disposition'] = cd 
    response.mimetype='text'
    return response
    
 
 
@app.route('/translate', methods = ['POST','GET'])
def translate():
    global data_frame
    if request.method == 'POST':
        if request.form.get('action1') == 'Download translated sequences':
            return redirect(url_for('download_translated_sequences'))
        if request.form.get('action2') == 'Analyze and Annotate':
            return redirect(url_for('name'))
        if request.form.get('action3') == 'Validate':
            return redirect(url_for('validate_results'))
    n = session.get('user_name',None)
    index =  list(range(1,len(sequences)+1)) 
    if(n[0]=='W'):
        protein_sequences = wipo_file[wipo_file['patent_number']== n]['Protein'].to_list()
    else:
        protein_sequences =uspto_file[uspto_file['patent_number']==n]['Protein'].to_list()
    data_frame = pd.DataFrame(list(zip(index, sequences,protein_sequences)),columns =['Index', 'Sequence','Translated Sequence'])
    return render_template('translated_sequences.html', column_names = data_frame.columns.values, row_data = list(data_frame.values.tolist()),zip = zip)
   
# --------- Prachi's edit------------#
@app.route('/download_translated_sequences',methods = ['POST','GET'])
def download_translated_sequences():
    response = make_response(data_frame.to_csv())
    cd = 'attachment; filename=translated_sequences.csv'
    response.headers['Content-Disposition'] = cd 
    response.mimetype='text/csv'
    return response
    
@app.route('/name',methods = ['POST','GET'])
def name():
    global seq_list
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
 
    return render_template('item.html',list_new = seq_list, user = n)
 
 
@app.route('/validate_results', methods = ['POST','GET'])
def validate_results():
    if request.method == 'POST':
        validate_sequence = request.form['seq']
        session['new_seq'] = validate_sequence
        if request.form.get('action1') == 'Validate for Processed Data':
            return redirect(url_for('actual_validation'))
        if request.form.get('action2') == 'Validate for Unprocessed Data':
            return redirect(url_for('unprocessed_validation'))
    n = session.get('user_name',None)
    if(n[0]=='W'):
        sequences_df = wipo_file[['patent_number','Sequence']]
        sequences_df = sequences_df[sequences_df['patent_number'] == n] 
        sequences = sequences_df['Sequence'].to_list()
    else:
        sequences_df = uspto_file[['patent_number','Sequence']]
        sequences_df = sequences_df[sequences_df['patent_number'] == n] 
        sequences = sequences_df['Sequence'].to_list()
    return render_template('validation.html',list_new = sequences, user = n)
 
@app.route('/actual_validation',methods = ['POST','GET'])
def actual_validation():
    validate_sequence = session.get('new_seq',None)
    if validate_sequence in sequences:
        statement = "The given sequence is present in the database on number " + str(sequences.index(validate_sequence)+1)
        return statement
    else:
        return "The given sequence is not present in the database for selected patent id."
    #seq_one = session.get('seq1',None)
    #alignment(seq_one,validate_sequence)
    #path = "/mnt/alignment.txt"
    #return send_file(path, as_attachment = True)
 
 
    #return render_template('validation_results.html', seq_one = validate_sequence, score = score, var = alignment)
 
@app.route('/unprocessed_validation', methods = ['POST','GET'])
def unprocessed_validation():
    validate_sequence = session.get('new_seq',None)
    if (validate_sequence[0].isupper()):
        result = process_protein(validate_sequence)
    else:
        result = process_string(validate_sequence)
    if result in sequences:
        statement = "The given sequence is present in the database on number " + str(sequences.index(result)+1)
        return statement
    else:
        return "The sequence is not present in the database for selected patent id."
    
    
@app.route('/display', methods = ['POST','GET'])
def display():
    global df
    if request.method == 'POST':
        return redirect(url_for('download_sequence'))
    it = session.get('chosen',None)
    print(it)
    # ---------- Prachi's edit ------------ #
    if(it == 'choose all'):
        print(seq_list)
        fasta_file_string = ",".join(seq_list)
        print(fasta_file_string)
        # Replacing the string with new line and > character for FASTA format #
        fasta_file = fasta_file_string.replace(",","\n>\n")
        print(fasta_file)
        # Adding > for the first sequence
        fasta_file = ">\n" + fasta_file
        i = 1
        fasta = ""
        for line in fasta_file:
            if(line == ">"):
                # To give numbers to each sequence #
                new_line = line + str(i)
                fasta += new_line
                i += 1
            else:
                fasta += line
        print(fasta)
        sequence = [{"name":"test","text": fasta}]
        bid = run_vel_blast(sequence,token,algo="blastp",collection="IC_TOXIN-All",evalue=10,nhits=250,word_size=3)
        print(bid)
        res=get_result(bid,token,30)
        print(res)
        with open("blas_res6.csv","w") as file:
            file.write(res.text)
        res_file=pd.read_csv("blas_res6.csv", header=None) 
        res_file.to_csv('blas_res6.csv', header=['qacc', 'sacc', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'], index=False)
        df= pd.read_csv("blas_res6.csv") 
        os.remove('blas_res6.csv')
    else:
        sequence=[{"name":"test","text":it}]
        print(sequence)
        bid=run_vel_blast(sequence,token,algo="blastp",collection="IC_TOXIN-All",
        evalue=10, nhits=250, word_size=3)
        print(bid)
        res=get_result(bid,token,30)
        print(res)
        with open("blas_res5.csv","w") as file:
            file.write(res.text)
        res_file=pd.read_csv("blas_res5.csv", header=None) 
        res_file.to_csv('blas_res5.csv', header=['qacc', 'sacc', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'], index=False)
        df= pd.read_csv("blas_res5.csv") 
        os.remove('blas_res5.csv')
    return render_template('results.html',seq = it,column_names = df.columns.values,row_data = list(df.values.tolist()),zip = zip)
 
 
# ------Prachi's edit ------- #
@app.route("/download_sequence", methods= ['GET','POST'])
def download_sequence():
    response = make_response(df.to_csv())
    cd = 'attachment; filename=mycsv.csv'
    response.headers['Content-Disposition'] = cd 
    response.mimetype='text/csv'
 
    return response
    
 

if __name__ == '__main__':
	app.run()
