from flask import Flask
from flask import render_template, url_for,redirect
from flask import request
from flask import session

app = Flask(__name__)
app.secret_key = "hello"


@app.route('/',methods = ['POST','GET'])
def index():
	list = ['Drinks','Vegetables','Fruits']
	if request.method == 'POST':
		n = request.form['n1']
		session['user_name'] = n
		#return user_name
		return redirect(url_for('name'))	
	return render_template('list.html', list = list)

@app.route('/name',methods = ['POST','GET'])
def name():
	if request.method == 'POST':
		it = request.form['n2']
		session['chosen'] = it
		return  redirect(url_for('display'))
	n = session.get('user_name',None)
	if n == 'Fruits':
		food_items = ['apple','orange','banana']
		return render_template('item.html',list_new = food_items,user = n)
	elif n == 'Vegetables':
		food_items = ['Okra','Potato','Capsicum']
		return render_template('item.html',list_new = food_items, user = n)
	elif n == 'Drinks':
		food_items = ['Orange juice','Lemon juice','Strawberry shake']
		return render_template('item.html', list_new = food_items, user = n)
	else:
		return 'OOPS'

@app.route('/display')
def display():
	return 'Great'

if __name__ == '__main__':
	app.run()
