#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<cmath>
#include<windows.h>
typedef double* douptr;
using namespace std;
int err = 0, ans_num = 1;
class matrix {
public:
	matrix(int row, int col): m(row), n(col) {
		mat = new douptr[m];
		for (int i = 0; i < m; i++) {
			mat[i] = new double[n];
		}
	}
	matrix():m(0),n(0),mat(nullptr){
	}
	matrix(const matrix& A) :m(A.m),n(A.n),name(A.name){
		mat = new douptr[m];
		for (int i = 0; i < m; i++) {
			mat[i] = new double[n];
			for (int j = 0; j < n; j++) {
				mat[i][j] = A.mat[i][j];
			}
		}

	}
	~matrix() {
		for (int i = 0; i < m; i++) {
			delete[] mat[i];
		}
		delete[] mat;
	}
	matrix& operator =(const matrix& R) {
		if (this == &R) {
			return *this;
		}
		for (int i = 0; i < m; i++) {
			delete[] mat[i];
		}
		delete[] mat;
		m = R.m;
		n = R.n;
		name = R.name;
		mat = new douptr[m];
		for (int i = 0; i < m; i++) {
			mat[i] = new double[n];
			for (int j = 0; j < n; j++) {
				mat[i][j] = R.mat[i][j];
			}
		}
		return *this;
	}
	void label(string A) {
		name = A;
	}
	friend ostream& operator<<(ostream& out, matrix& M) {
		for (int i = 0; i < M.m; i++) {
			//answer label
			if (i == M.m / 2) {
				out << M.name << " = ";
			}
			out << '\t';
			//lu corner
			if (i == 0) {
				out << "┌  ";
			}
			//ld corner
			else if (i == M.m - 1 && M.m != 1) {
				out << "└  ";
			}
			//l side
			else {
				out << "│  ";
			}

			for (int j = 0; j < M.n; j++) {
				cout << M.mat[i][j];
				if (j != M.n - 1) {
					cout << '\t';
				}
			}
			//ru corner
			if (i == 0) {
				out << "  ┐";
			}
			else if (i == M.m - 1 && M.m != 1) {
				out << "  ┘";
			}
			else {
				out << "  │";
			}
			if (i == M.m - 1) {
				cout << '(' << M.m << "×" << M.n << ')';
			}
			cout << endl;
		}
		return out;
	}
	friend istream& operator>>(istream& in, matrix& M){
		for (int i = 0; i < M.m; i++) {
			for (int j = 0; j < M.n; j++) {
				in >> M.mat[i][j];
			}
		}
		return in;
	}
	friend matrix operator +(matrix A, matrix B){
		matrix temp;
		if (A.m != B.m || A.n != B.n) {
			throw ("错误：矩阵"+ A.name+ "与"+ B.name+ "不同型，无法相加");
		}
		temp = A;
		for (int i = 0; i < A.m; i++) {
			for (int j = 0; j < A.n; j++) {
				temp.mat[i][j] += B.mat[i][j];
			}
		}
		return temp;
	}
	friend matrix operator -(matrix A, matrix B){
		matrix temp;
		if (A.m != B.m || A.n != B.n) {
			throw ("错误：矩阵" + A.name + "与" + B.name + "不同型，无法相加");
		}
		temp = A;
		for (int i = 0; i < A.m; i++) {
			for (int j = 0; j < A.n; j++) {
				temp.mat[i][j] -= B.mat[i][j];
			}
		}
		return temp;
	}
	friend matrix operator *(matrix A, matrix B){
		matrix temp(A.m, B.n);
		if (A.n != B.m) {
			throw ( "错误：存在无法相乘的矩阵" + A.name + "和" + B.name );
		}
		double sum = 0.0;
		for (int i = 0; i < A.m; i++) {
			for (int j = 0; j < B.n; j++) {
				for (int k = 0; k < A.n; k++) {
					sum += A.mat[i][k] * B.mat[k][j];
				}
				temp.mat[i][j] = sum;
				sum = 0;
			}
		}
		return temp;
	}
	friend matrix operator *(double B, matrix A){
		matrix temp = A;
		for (int i = 0; i < A.m; i++) {
			for (int j = 0; j < A.n; j++) {
				temp.mat[i][j] *= B;
			}
		}
		return temp;
	}
	friend matrix operator *(matrix A, double B){
		matrix temp = A;
		for (int i = 0; i < A.m; i++) {
			for (int j = 0; j < A.n; j++) {
				temp.mat[i][j] *= B;
			}
		}
		return temp;
	}
	friend matrix transpose(matrix A){
		matrix temp(A.n, A.m);
		for (int i = 0; i < A.m; i++) {
			for (int j = 0; j < A.n; j++) {
				temp.mat[j][i] = A.mat[i][j];
			}
		}
		return temp;
	}
	friend matrix cofactor(matrix A, int r, int c) {
		int size = A.m - 1;
		if (size <= 0) {
			throw "矩阵太小，无余子式";
		}
		matrix temp(size, size);
		for (int i = 0; i < size + 1; i++) {
			for (int j = 0; j < size + 1; j++) {
				if (i > r && j < c) {
					temp.mat[i-1][j]= A.mat[i][j];
				}
				else if (i<r && j>c) {
					temp.mat[i][j-1] = A.mat[i][j];
				}
				else if (i>r&&j>c) {
					temp.mat[i-1][j-1] = A.mat[i][j];
				}
				else if(i<r&&j<c){
					temp.mat[i][j] = A.mat[i][j];
				}
			}
		}
		return temp;
	}
	friend double determinant(matrix A) {
		double result = 0;
		if (A.m != A.n) {
			throw "非方阵，无法求行列式！";
		}
		if (A.m == 1) {
			result = A.mat[0][0];
			return result;
		}
		if (A.m == 2) {
			result = A.mat[0][0] * A.mat[1][1] - A.mat[0][1] * A.mat[1][0];
			return result;
		}
		int one = 1;
		for (int i = 0; i < A.m; i++) {
			result += one * A.mat[0][i] * determinant(cofactor(A, 0, i));
			one = -one;
		}
		return result;
	}

	string nameit(){
		return name;
	}
	void setij(int i, int j, double value = 0.0){
		mat[i][j] = value;
	}
private:
	douptr* mat;
	int m, n;
	string name;
};
vector<matrix> S;
void create(int m, int n, string name) {
	matrix temp(m, n);
	cin >> temp;
	S.push_back(temp);
	S[S.size()-1].label(name);
}
int str_to_int(string str) {
	for (int i = 0; i < str.length(); i++) {
		if (str[i] < '0' || str[i]>'9') {
			throw "矩阵大小应填格式正确的整数";
		}
	}
	int result = 0;
	for (int i = str.length() - 1; i >= 0; i--) {
		result += static_cast<int>(str[i] - '0') * pow(10, str.length() - 1 - i);
	}
	return result;
}
double str_to_dou(string str) {
	for (int i = 0; i < str.length(); i++) {
		if (str[i] < '0' || str[i]>'9'||str[i]=='.') {
			throw "存在形式上错误的数字";
		}
	}
	bool point = false;
	int count = 0;
	double result=0.0;
	int len = str.length();
	for (int i = 0; i < len; i++) {
		if (point) {
			count++;
		}
		if (str[i] == '.') {
			point = true;
			str.erase(i, 1);
		}
	}
	result = static_cast<double>(str_to_int(str)) / pow(10, count);
	return result;
}
string timename() {
	SYSTEMTIME sys;
	GetLocalTime(&sys);
	string name = to_string(sys.wYear) +"-" + to_string(sys.wMonth) +"-"+ to_string(sys.wDay) + "-" + to_string(sys.wHour) + "-" + to_string(sys.wMinute) + "-" + to_string(sys.wSecond);
	return name;
}
void savefile(string name ) {
	ofstream file(name + ".txt");
	if (!file) {
		throw "文件创建失败";
	}
	for (int i = 0; i < S.size(); i++) {
		file << S[i].nameit() <<"：" << endl;
		file << S[i];
	}
	file.close();
	cout << "文件已保存为" << name << ".txt" << endl;
}
void rules() {
	cout << "语法：\n定义一个m×n矩阵：define 矩阵名(m,n)" << endl
		<< "特别地，定义n×n单位矩阵：define unitm 矩阵名(n)" << endl
		<< "定义n×n的零矩阵：define zerom 矩阵名(n)" << endl
		<< "显示一个矩阵：show 矩阵名" << endl
		<< "计算：calculate 表达式" << endl
		<< "转置：transpose 矩阵名，或define 矩阵名=transposed 已有的某矩阵。" << endl
		<< "注：表达式仅支持+、-（不能在最前）"
		<< "*或·（乘），其中乘仅限数字与矩阵和矩阵与矩阵。" << endl
		<< "tips：向量即是行或列为1的矩阵。" << endl
		<< "求行列式：det 矩阵名" << endl
		<< "保存文件（自动命名）：save " << endl
		<< "自定义文件名：save as 文件名";
}
void commander() {
	string cmd;
	getline(cin, cmd);
	if (cmd == "rules") {
		rules();
		return;
	}
	if (cmd.substr(0, 6) == "define") {
		string name;
		int m, n;
		if (cmd.substr(7, 5) == "unitm") {
			if (cmd.find('(', 13) == string::npos) {
				throw "应有\"(\"为界";
			}
			if (cmd.find(')', 13) == string::npos) {
				throw "应有\")\"为界";
			}
			name = cmd.substr(13, cmd.find('(', 13) - 13);
			m = str_to_int(cmd.substr(cmd.find('(', 13) + 1, cmd.find(')', 13) - cmd.find('(', 13) - 1));
			matrix temp(m, m);
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < m; j++) {
					if (i == j) {
						temp.setij(i, j, 1.0);
					}
					else {
						temp.setij(i, j);
					}
				}
			}
			S.push_back(temp);
			S[S.size()-1].label(name);
			cout << "已完成定义" << endl;
			return;
		}
		if (cmd.substr(7, 5) == "zerom") {
			if (cmd.find('(', 13) == string::npos) {
				throw "应有\"(\"为界";
			}
			if (cmd.find(')', 13) == string::npos) {
				throw "应有\")\"为界";
			}
			name = cmd.substr(13, cmd.find('(', 13) - 13);
			m = str_to_int(cmd.substr(cmd.find('(', 13) + 1, cmd.find(')', 13) - cmd.find('(', 13) - 1));
			matrix temp(m, m);
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < m; j++) {
					temp.setij(i, j);
				}
			}
			S.push_back(temp);
			S[S.size()-1].label(name);
			cout << "已完成定义" << endl;
			return;
		}
		//below has a problem...
		if (cmd.find("=transposed")!=string::npos) {
			//cmd.substr(cmd.find("="), 11) == "=transposed"
			string toseek = cmd.substr(cmd.find("=")+12);
			name = cmd.substr(7, cmd.find("=") - 7);
			bool find = false;
			int i = 0;
			for ( i= 0; i < S.size(); i++) {
				if (S[i].nameit() == toseek) {
					find = true;
					break;
				}
			}
			if (find) {
				S.push_back(transpose(S[i]));
				S[S.size()-1].label(name);
				cout << "定义成功。" << endl;
				return;
			}
			else {
				throw ("未找到矩阵" + toseek);
			}
		}
		name = cmd.substr(7, cmd.find('(', 7) - 7);
		m = str_to_int(cmd.substr(cmd.find('(', 7) + 1, cmd.find(',', cmd.find('(', 7)) - cmd.find('(', 7) - 1));
		n = str_to_int(cmd.substr(cmd.find(',', cmd.find('(', 7)) + 1, cmd.find(')', cmd.find(',', cmd.find('(', 7))) - cmd.find(',', cmd.find('(', 7)) - 1));
		create(m, n, name);
	}
	if (cmd.substr(0, 4) == "show") {
		bool find = false;
		for (int i = 0; i < S.size(); i++) {
			if (S[i].nameit() == cmd.substr(5)) {
				cout << S[i];
				find = true;
				return;
			}
		}
		if (!find) {
			throw ("未找到矩阵" + cmd.substr(5));
		}
	}
	if (cmd.substr(0, 9) == "calculate") {
		vector<string> names;
		vector<double> nums;
		vector<int> mat_index;
		vector<int> pm_index;
		cmd.erase(0, 10);
		for (int i = 0; i < cmd.length(); i++) {
			if (cmd[0] == '-') {
				throw "不能以“-”开头，建议定义合适的零矩阵并使用O-A";
			}
		}
		int i, j, start = 0;
		string Temp;
		for (i = 0; i < cmd.length(); i++) {
			if (cmd[i] == '+' || cmd[i] == '-' || cmd[i] == '*' || cmd[i] == '·') {
				names.push_back(Temp);
				if (cmd[i] == '+') {
					names.push_back("+");
				}
				else if (cmd[i] == '-') {
					names.push_back("-");
				}
				else {
					names.push_back("*");
				}
				Temp.clear();
			}
			else {
				Temp += cmd[i];
			}
		}

		names.push_back(Temp);
		//to analyze the composition:
		for (i = 0; i < names.size(); i++) {
			bool is_mat = false;
			if (names[i] != "+" && names[i] != "-" && names[i] != "*") {
				//to see whether it's a matrix:
				for (j = 0; j < S.size(); j++) {
					if (S[j].nameit() == names[i]) {
						mat_index.push_back(j);
						is_mat = true;
						names[i] = "a matrix";
						break;
					}
				}
				//to figure out whether it's a number:
				if (!is_mat) {
					nums.push_back(str_to_dou(names[i]));
					names[i] = "a number";
				}
			}
			else {
				if (names[i] != "*") {
					pm_index.push_back(i);
				}
			}
		}
		//to complete the calculation:
		int num_used = 0, mat_used = 0;
		bool is_initialized = false;
		double num_temp = 1.0;
		pm_index.push_back(names.size());
		vector<matrix> temp(pm_index.size());
		for (i = 0; i < pm_index.size(); i++) {
			for (j = start; j < pm_index[i]; j++) {
				if (names[j] == "a number") {
					num_temp *= nums[num_used++];
				}
				if (names[j] == "a matrix") {
					if (is_initialized) {
						temp[i] = temp[i] * S[mat_index[mat_used]];
					}
					else {
						temp[i] = S[mat_index[mat_used++]];
						is_initialized = true;
					}
				}
			}
			temp[i] = temp[i] * num_temp;
			num_temp = 1.0;
			is_initialized = false;
			start = pm_index[i] + 1;
		}
		matrix Ans = temp[0];
		for (i = 0; i < pm_index.size() - 1; i++) {
			if (names[pm_index[i]] == "+") {
				Ans = Ans + temp[i + 1];
			}
			else {
				Ans = Ans - temp[i + 1];
			}
		}

		S.push_back(Ans);
		S[S.size() - 1].label("Ans" + to_string(ans_num++));
		cout << S[S.size() - 1];
		cout << "此答案已保存为" << S[S.size() - 1].nameit() << endl;
		return;
	}
	if (cmd.substr(0, 9) == "transpose") {
		cmd.erase(0, 10);
		bool find = false;
		int i;
		for (i = 0; i < S.size(); i++) {
			if (S[i].nameit() == cmd) {
				find = true;
				S[i] = transpose(S[i]);
				break;
			}
		}
		if (find) {
			cout << S[i].nameit() << "已完成转置。" << endl;
			return;
		}
		else {
			throw ("未找到矩阵" + cmd);

		}
		if (cmd.substr(0, 3) == "det") {
			return;
		}
	}
	if (cmd.substr(0, 3) == "det") {
		string name = cmd.substr(4);
		bool find = false;
		for (int i = 0; i < S.size(); i++) {
			if (S[i].nameit() == name) {
				find = true;
				cout << "= " << determinant(S[i]) << endl;
				return;
			}
		}
		if (!find) {
			throw ("未找到矩阵" + name);
		}
		return;
	}
	if (cmd.substr(0, 4) == "save") {
		if (cmd == "save") {
			savefile(timename());
			return;
		}
		if (cmd.substr(5) == "as") {
			cout << "请输入文件名：";
			string name;
			getline(cin, name);
			savefile(name);
			return;
		}
		bool find = false;
		for (int i = 0; i < S.size(); i++) {
			if (S[i].nameit() == cmd.substr(5)) {
				savefile(S[i].nameit());
				find = true;
				return;
			}
		}
		if (!find) {
			throw ("未找到矩阵" + cmd.substr(5));
		}
	}
		if (cmd == "save") {
			savefile(timename());
			return;
		}
		if (cmd.substr(5,2) == "as") {
			cout << "请输入文件名(不用加.txt)：";
			string name=cmd.substr(8);
			savefile(name);
			return;
		}
}
int main() {
	cout << "【TIPS】输入rules以了解使用方法。" << endl;
	while (true) {
		try {
			commander();
		}
		catch (string Errormessage) {
			cout << Errormessage << endl;
		}	
	 }
	return 0;
}