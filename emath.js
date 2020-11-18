function ematrix(rows,cols){
    this.rows = rows;
    this.cols = cols;
    this.data = [];
    for(var i = 0;i<rows;i++)
	this.data.push(new Array(cols));

    this.zeros = function(){
	for(var i = 0;i<this.rows;i++)
	    for(var j = 0;j<this.cols;j++)
		this.data[i][j] = 0;
    }
    this.ones = function(){
	for(var i = 0;i<this.rows;i++)
	    for(var j = 0;j<this.cols;j++)
		this.data[i][j] = 1;
    }
    this.random = function(){	
	for(var i = 0;i<this.rows;i++)
	    for(var j = 0;j<this.cols;j++)
		this.data[i][j] = Math.random();
    }
    this.eye = function(){
	for(var i = 0;i<this.rows;i++){
	    for(var j = 0;j<this.cols;j++){
		if(i==j)
		    this.data[i][j] = 1;
		else
		    this.data[i][j] = 0;
	    }
	}
    }
    this.to_string = function(){
	var s = "";
	for(var i = 0;i<this.rows;i++){
	    for(var j = 0;j<this.cols;j++)
		s += this.data[i][j] + " ";
	    s += "\n";
	}
	return s;
    }
    this.add = function(M){
	var out = new ematrix(this.rows,cols);
	for(var i = 0;i<this.rows;i++){
	    for(var j = 0;j<this.cols;j++)
		out.data[i][j] = this.data[i][j] + M.data[i][j];
	}
	return out;
    }
    this.diff = function(M){
	var out = new ematrix(this.rows,cols);
	for(var i = 0;i<this.rows;i++){
	    for(var j = 0;j<this.cols;j++)
		out.data[i][j] = this.data[i][j] - M.data[i][j];
	}
	return out;
    }
    this.mul = function(M){
	var out = new ematrix(this.rows,cols);
	for(var i = 0;i<this.rows;i++){
	    for(var j = 0;j<this.cols;j++)
		out.data[i][j] = this.data[i][j] * M.data[i][j];
	}
	return out;
    }
    this.div = function(M){
	var out = new ematrix(this.rows,cols);
	for(var i = 0;i<this.rows;i++){
	    for(var j = 0;j<this.cols;j++)
		out.data[i][j] = this.data[i][j] / M.data[i][j];
	}
	return out;
    }
    this.dot = function(M){
	var out = new ematrix(this.rows,M.cols);
	var s;
	for(var i = 0;i<this.rows;i++){
	    for(var j = 0;j<M.cols;j++){
		s = 0;
		for(var k = 0;k<this.cols;k++)
		    s += this.data[i][k] * M.data[k][j];
		out.data[i][j] = s;
	    }
	}
	return out;
    }
    this.scale = function(num){
	var out = new ematrix(this.rows,cols);
	for(var i = 0;i<this.rows;i++){
	    for(var j = 0;j<this.cols;j++)
		out.data[i][j] = this.data[i][j] * num;
	}
	return out;
    }
    this.trace = function(){
	var s = 0;
	for(var i = 0;i<this.rows;i++)
	    s+=this.data[i][i];
	return s;
    }
    this.diag = function(){
	var out = new ematrix(this.rows,this.cols);
	out.zeros();
	for(var i = 0;i<this.rows;i++)
	    out.data[i][i] = this.data[i][i];
	return out;
    }
    this.transpose = function(){
	var out = new ematrix(this.cols,this.rows);
	for(var i = 0;i<this.rows;i++)
	    for(var j = 0;j<this.cols;j++)
		out.data[j][i] = this.data[i][j];
	return out;
    }
    this.lu = function(){
	var L = new ematrix(this.rows,this.cols);
	var U = new ematrix(this.rows,this.cols);
	var s;
	L.eye();
	U.zeros();
	for(var i = 0;i<this.rows;i++){
	    for(var j = i;j<this.cols;j++){
		s = 0;
		for(var k = 0;k<i;k++)
		    s += L.data[i][k]*U.data[k][j];
		U.data[i][j] = this.data[i][j] - s;
	    }
	    for(var j = i+1;j<this.cols;j++){
		s = 0;
		for(var k = 0;k<i;k++)
		    s += L.data[j][k]*U.data[k][i];
		L.data[j][i] = (this.data[j][i] - s)/U.data[i][i];
	    }
	}
	return [L,U];
    }
    this.qr = function(){
	var Q = new ematrix(this.rows,this.cols);
	var R = new ematrix(this.rows,this.cols);
	var tmp = new ematrix(this.rows,1);
	var s;
	R.zeros();
	Q.zeros();
	for(var i = 0;i<this.rows;i++){
	    tmp.zeros();
	    for(var j = 0;j<i;j++){
		s = 0;
		for(var k = 0;k<this.rows;k++)
		    s += this.data[k][i] * Q.data[k][j];
		R.data[j][i] = s;
		for(var k = 0;k<this.rows;k++)
		    tmp.data[k][0] += R.data[j][i]*Q.data[k][j];
	    }
	    s = 0;
	    for(var j = 0;j<this.rows;j++){
		Q.data[j][i] = this.data[j][i] - tmp.data[j][0];
		s += Q.data[j][i]*Q.data[j][i];
	    }
	    s = Math.sqrt(s);
	    R.data[i][i] = s;
	    for(var j = 0;j<this.rows;j++)
		Q.data[j][i] /= s;
	}
	return [Q,R];
    }
    this.det = function(){
	var lu;
	var s = 1;
	lu = this.lu();
	for(var i = 0;i<this.rows;i++)
	    s *= lu[1].data[i][i];
	return s;
    }
    this.inv = function(){
	var out = new ematrix(this.rows,this.cols);
	var z = new ematrix(this.rows,1);
	var lu = this.lu();
	var L,U;
	var s;
	L = lu[0];
	U = lu[1];
	for(var j = 0;j<this.cols;j++){
	    for(var i = 0;i<this.rows;i++){
		s = 0;
		for(var k = 0;k<i;k++)
		    s += L.data[i][k] * z.data[k][0];
		if(i==j)
		    z.data[i][0] = 1 - s;
		else
		    z.data[i][0] = - s;
	    }
	    for(var i = this.rows - 1;i>=0;i--){
		s = 0;
		for(var k = i+1;k<this.rows;k++)
		    s += U.data[i][k] * out.data[k][j];
		out.data[i][j] = (z.data[i][0] - s)/U.data[i][i];
	    }
	}
	return out;
    }
    this.eig = function(maxIt){
	var qr;
	var eig_val,eig_vec;
	eig_val = this;
	eig_vec = new ematrix(this.rows,this.cols);
	eig_vec.eye();
	for(var i = 0;i<maxIt;i++){
	    qr = eig_val.qr();
	    eig_vec = eig_vec.dot(qr[0]);
	    eig_val = qr[1].dot(qr[0]);
	}
	return [eig_val,eig_vec];
    }
}

function epoly(coeffs){
    this.degree = coeffs.length - 1
    this.coeffs = coeffs;

    this.to_string = function(){
	var s = "";
	s += this.coeffs[0];
	for(var i =1;i<=this.degree;i++)
	    s += "+"+this.coeffs[i] + "x^"+ i;
	return s;
    }
    this.eval = function(x){
	var x0 = 1;
	var s = 0;
	for(var i = 0;i<=this.degree;i++){
	    s += this.coeffs[i]*x0;
	    x0 *= x;
	}
	return s;
    }
    this.derivate = function(){
	var out = [];
	for(var i = 1;i<=this.degree;i++)
	    out.push(i*this.coeffs[i]);
	return new epoly(out);
    }
    this.integrate = function(){
	var out = [];
	out.push(0);
	for(var i = 0;i<=this.degree;i++)
	    out.push(this.coeffs[i]/(i+1));
	return new epoly(out);
    }
    this.add = function(p){
	var out = [];
	var s;
	var d = this.degree>p.degree?this.degree:p.degree;
	for(var i = 0;i<=d;i++){
	    s = 0;
	    if(i<=this.degree)
		s += this.coeffs[i];
	    if(i<=p.degree)
		s += p.coeffs[i];
	    out.push(s);
	}
	return new epoly(out);
    }
    this.diff = function(p){
	var out = [];
	var s;
	var d = this.degree>p.degree?this.degree:p.degree;
	for(var i = 0;i<=d;i++){
	    s = 0;
	    if(i<=this.degree)
		s += this.coeffs[i];
	    if(i<=p.degree)
		s -= p.coeffs[i];
	    out.push(s);
	}
	return new epoly(out);
    }
    this.mul = function(p){
	var d = this.degree + p.degree + 1;
	var out = new Array(d);
	for(var i = 0;i<d;i++)
	    out[i] = 0;
	for(var i = 0;i<=this.degree;i++){
	    for(var j = 0;j<=p.degree;j++)
		out[i+j] += this.coeffs[i]*p.coeffs[j];
	}
	return new epoly(out);
    }
    this.div = function(p){
	var q = [];
	var r = [];
    }
    this.scale = function(num){
	var out = [];
	for(var i = 0;i<=this.degree;i++)
	    out.push(this.coeffs[i]*num);
	return new epoly(out);
    }
}

let a = new epoly([1,2,1]);
let b = new epoly([1,1]);
console.log(a.to_string());
console.log(a.diff(b).to_string());
