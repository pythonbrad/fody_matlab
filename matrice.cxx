#include <iostream>
#include <string>
#define MAX_ORDER 100
#define MAX_NB_MATRIX 50
#define ERROR_MATRIX_NOT_FOUND(tag) cout << "Error: matrix " << data << " unknowed;" << endl
#define MAX_PRECISION 1e-15
// return 0 if value in [-MAX_PRECISION, MAX_PRECISION]
#define FIX_PRECISION(varname) varname = (varname < MAX_PRECISION || -varname > -MAX_PRECISION) ? 0 : varname
using namespace std;

class Matrix {
    double data[MAX_ORDER*MAX_ORDER];
    public:
		string tag="";
        unsigned short order;
        Matrix() {
            this->order = 0;
        };
        Matrix(unsigned char x) {
            this->order = x;
        };
        Matrix operator + (Matrix matrix) {
            unsigned short row, col;
            Matrix temp(this->order);
            for (row = 0; row < order; ++row) {
                for (col = 0; col < order; ++col) {
                    temp.write(row, col, this->read(row, col) + matrix.read(row, col));
                };
            };
            return temp;
        };
        Matrix operator * (Matrix matrix) {
            unsigned short row, col, c = 0;
            Matrix temp(this->order);
            for (row = 0; row < order; ++row) {
                for (col = 0; col < order; ++col) {
                    temp.write(row, col, 0.0);
                    for (c = 0; c < order; ++c) {
                        temp.write(row, col, temp.read(row, col) + this->read(row, c) * matrix.read(c, col));
                    };
                };
            };
            return temp;
        };
        bool operator == (Matrix matrix) {
			unsigned short row, col;
			double a,b;
			for (row = 0; row < order; ++row) {
                for (col = 0; col < order; ++col) {
					a = this->read(row,col);
					b = matrix.read(row, col);
					// This fix the problem of precision
					FIX_PRECISION(a);
					FIX_PRECISION(b);
					if (a != b) {
						return 0;
					};
				};
			};
			return 1;
		};
        Matrix operator - (Matrix matrix) {
            return operator + (matrix * -1);
        };
        Matrix operator * (double k) {
            Matrix temp(this->order);
            temp.build_from_scalar(k);
            return operator * (temp);
        };
        Matrix operator / (Matrix matrix) {
            return *this * matrix.invert();
        };
        Matrix operator / (double k) {
            Matrix temp(this->order);
            temp.build_from_scalar(k);
            return operator / (temp);
        };
        void write(unsigned short row, unsigned short col, double value) {
            this->data[row*(this->order)+col] = value;
        };
        double read(unsigned short row, unsigned short col) {
            return this->data[row*(this->order)+col];
        };
        void build_from_scalar(double scalar) {
            unsigned short row, col;
            for (row = 0; row < this->order; ++row) {
                for (col = 0; col < this->order; ++col) {
                    this->write(col, row, (row==col) ? scalar : 0);
                };
            };
        };
        void build_from_input() {
            unsigned short row;
            unsigned short col;
            double value;
            string entry;
            for (row=0; row < order; row++) {
                for (col=0; col < order; col++) {
                    printf("Enter the value (%i,%i): ", row, col);
                    cin >> entry;
                    value = this->get_decimal(entry);
                    // iu verify if error
                    if (value == 0 && entry != "0") {
						cout << "Enter R (retry) or another key to cancel: ";
						cin >> entry;
						if (entry == "R") {
							col--;
							continue;
						} else break;
					} else {
						this->write(row, col, value);
					}
                };
            };
        };
        bool build_from_string(string code) {
            unsigned short i = 1;
            unsigned short next_id = 0;
            string token;
            while (1) {
                if (code[i] == ',' || code[i] == '}') {
                    this->data[next_id] = get_decimal(token);
                    next_id++;
                    token = "";
                    if (code[i] == '}') break;
                } else {
                    token += code[i];
                };
                i++;
            };
            for (i=0; i*i!=next_id; i++) {
				if (i>next_id+1) {
					//cout << "The matrix data not correspond to a square matrix";
					this->build_from_scalar(0);// I clean the matrix
					return 0;
				};
			};
			if (this->order!=i) {
				//cout << "The nb of elements should be " << (this->order)*(this->order);
				this->build_from_scalar(0);// I clean the matrix
				return 0;
			};
			return 1;
        };
        double trace() {
            double result = 0.0;
            unsigned short i;
            for (i = 0; i < this->order; ++i) {
                result += this->read(i, i);
            };
            return result;
        };
        double determinant() {
            unsigned short row, col;
            double result, sum_of_product, diff_of_product;
            if (this->order > 2) {
                result = 0.0;
                for (row = 0; row < this->order; ++row) {
                    sum_of_product = diff_of_product = 1;
                    for (col = 0; col < this->order; ++col) {
                        sum_of_product *= this->read((row+col)%this->order, col);
                        diff_of_product *= this->read((row+col)%this->order, this->order-1-col);
                    };
                    result += sum_of_product - diff_of_product;
                };
            } else if (order == 2) {
                result = this->read(0, 0)*this->read(1, 1) - this->read(1, 0)*this->read(0, 1);
            } else {
				return this->read(0,0);
			};
            return result;
        };
        Matrix cofactor() {
            unsigned short col, row;
            unsigned short col2, row2;
            Matrix Mresult(this->order);
            Matrix Mtemp(this->order-1);
            for (row = 0; row < this->order; ++row) {
                for (col = 0; col < this->order; ++col) {
                    for (row2 = 0; row2 < this->order; ++row2) {
						if (row2 == row) continue;
						for (col2 = 0; col2 < this->order; ++col2) {
							if(col2 == col) continue;
							// -1 because 1 row or col has been skipped
							Mtemp.write((row2 > row) ? (row2 - 1) : row2, (col2 > col) ? (col2 - 1) : col2, this->read(row2, col2));
						};
					};
					Mresult.write(row, col, ((((row+col) % 2) == 0) ? 1 : -1) * Mtemp.determinant());
                };
            };
            return Mresult;
        };
        Matrix transpose() {
            unsigned short row, col;
            Matrix temp(this->order);
            for (row = 0; row < this->order; ++row) {
                for (col = 0; col < this->order; ++col) {
                    temp.write(col, row, this->read(row, col));
                };
            };
            return temp;
        };
        Matrix adjacent() {
            return this->cofactor().transpose();
        };
        Matrix invert() {
            return this->adjacent() * (1/this->determinant());//cannot use "/" directly, because this operator use invert
        };
        void print() {
			// NB: IF PROBLEM OF PRINT, VERIFY THE TYPES, CONVERSION FIRTLY, double value (can be too small eg:1e-20)
            unsigned short row, col, c, ms=0, cs=0;
            double temp;
            // I get the max number of digits
            for (row = 0; row < this->order; ++row) {
                for (col = 0; col < this->order; ++col) {
                    temp = this->read(row, col);
                    cs = this->get_number_of_digits(temp);
                    if (cs > ms) ms = cs;
                };
            };
            for (row = 0; row < this->order; ++row) {
                for (col = 0; col < this->order; ++col) {
                    temp = this->read(row, col);
                    // I apply the max precision, to fix some errors
                    FIX_PRECISION(temp);
                    // I put the spaces in the unused digit spaces
                    // I don't add 1 for a negative number because the sign take it
                    for (c=0; c<ms-this->get_number_of_digits(temp)+1; ++c) {
						cout << " ";
					};
                    printf("%.2lf", temp);
                };
                // I go to the next line
                cout << endl;
            };
        };
        bool test() {
            Matrix matrixA(3), matrixB(3), matrixC(3), matrixD(3), matrixE(6), matrixF(6), matrixG(3), matrixH(4);
            if (matrixA.build_from_string("{-1,2,-3,2,1,0,4,-2,5}") && matrixB.build_from_string("{-5,4,-3,10,-7,6,8,-6,5}")) {
				matrixA.print();
				matrixA.invert().print();
				matrixC.build_from_scalar(1);
				matrixD.build_from_scalar(0);
				matrixE.build_from_scalar(-1.125);
				matrixF.build_from_scalar(1);
				matrixG.build_from_string("{-0.00345,2.345,-3.345,0,-0.0945,0.10345,4.345,-2.345,5.345}");
				matrixG.print();
				(matrixG/matrixG).print();
				matrixH.build_from_scalar(5);
				matrixH.adjacent().print();
				return (matrixA.invert() == matrixB) && (matrixA*matrixB == matrixC) && \
				(matrixA-matrixA == matrixD) && matrixE/matrixE == matrixF && \
				matrixA.get_decimal("-123.123") == -123.123 && matrixG/matrixG == matrixC && \
				get_number_of_digits(0.000001) == 1 && get_number_of_digits(-0.00000000001) == 2 && \
				get_number_of_digits(1) == 1 && get_number_of_digits(-1) == 2 && get_number_of_digits(-1.1) == 2;
			} else return 0;
        };
        unsigned short get_number_of_digits(double num) {
			unsigned short n = (num < 0) ? (num=-num, 1) : 0;
			n += (num < 1) ? 1 : 0;
			while (num >= 1) {
                num /= 10;
                n += 1;
            };
            return n;
		};
        double get_decimal(string text) {
            double decimal=0;
            double frac = 0;
            unsigned short i;
            bool s = (text[0] == '-') ? 1 : 0;
            for (i=s; text[i] != '\0'; ++i) {
                decimal *= 10;
                frac *= 10;
                switch (text[i]) {
                    case '0':
                        break;
                    case '1':
                        decimal += 1;
                        break;
                    case '2':
                        decimal += 2;
                        break;
                    case '3':
                        decimal += 3;
                        break;
                    case '4':
                        decimal += 4;
                        break;
                    case '5':
                        decimal += 5;
                        break;
                    case '6':
                        decimal += 6;
                        break;
                    case '7':
                        decimal += 7;
                        break;
                    case '8':
                        decimal += 8;
                        break;
                    case '9':
                        decimal += 9;
                        break;
                    case '.':
						frac = 1;
						decimal /= 10;
						break;
                    default:
                        cout << "error while the getting of decimal, char " << text[i] << " unknowed" << endl;
                        return 0;
                };
            };
            return decimal/(frac ? frac : 1)*(s ? -1.0 : 1.0);
        };
};

bool is_character(char c) {
    char chars[28] = "abcdefghijklmnopqrstuvwxyz$";
    int i;
    for (i = 0; i < 28; ++i) {
        if (chars[i] == c) {
            return 1;
        };
    };
    return 0;
};

bool is_digit(char c) {
	char digits[11] = "0123456789";
	int i;
    for (i = 0; i < 11; ++i) {
        if (digits[i] == c) {
            return 1;
        };
    };
    return 0;
};

bool is_space(char c) {
    return (c == ' ' || c == '\n' || c == '\t');
};

bool is_operator(char c) {
    char operators[7] = "+-/*&.";
    int i;
    for (i = 0; i < 7; ++i) {
        if (operators[i] == c) {
            return 1;
        };
    };
    return 0;
};

unsigned char is_special_operator(string data) {
	if (data == "trace") {
		return 'x';
	} else if (data == "adj") {
		return 'a';
	} else if (data == "cof") {
		return 'c';
	} else if (data == "inv") {
		return 'i';
	} else if (data == "trans") {
		return 't';
	} else if (data == "det") {
		return 'd';
	} else {
		return '\0';
	};
};

unsigned short getMatrixIdByTag (Matrix MATRIX[MAX_NB_MATRIX], string tag, unsigned short length) {
	unsigned short i;
	for (i=0; i<length; i++) {
		if (MATRIX[i].tag == tag) {
			return i;
		};
	};
	return MAX_NB_MATRIX;
}

int main() {
    Matrix MATRIX[MAX_NB_MATRIX];
    MATRIX[0].build_from_scalar(0);
    MATRIX[0].tag = '$';
	unsigned short next_matrix_id=1;//matrix counter
	unsigned short current_char=0;//string counter
    string cmd;
    char op;//operator
    unsigned short temp;
    double temp2;
    unsigned short MLS;//matrix to left side
    unsigned short state;
    // A:order
    // input
    // A:order = {...};
    // A + B;
    // $;
    // A.trace
    // 0 null
    // 10 initialization (order)
    // 11 initialization (value)
    // 12 initialization (data)
    // 13 saving
    // 20 operation
    // 1 error
    // 2 succeed
    string data;
    //TEST
    cout << (MATRIX[0].test() ? "Test OK" : "Test Failed") << endl;
    cout << "Welcome to FoDy MatLab (precision: " << MAX_PRECISION << ")" << endl;
    while (1) {
		cout << "> ";
		getline(cin, cmd);
		if (cmd != "exit") {
			// I clean data
			data = "";state = 0;
			for (current_char = 0; cmd[current_char] != '\0' && state != 1; current_char++) {
				if (cmd[current_char] == '#') {
					break;
				} else if (is_space(cmd[current_char])) {
					continue;
				} else if (cmd[current_char] == ';') {
					// Normally, I execute each instruction
					switch (state) {
						case 0:
							break;
						case 11:
						case 13:
							if (next_matrix_id < MAX_NB_MATRIX) {
								if (data == "$") {
									// Just copy
								} else if (data == "input") {
									MATRIX[0].build_from_input();
								} else if (state == 13) {
									if (!MATRIX[0].build_from_string('{'+data+'}')) {
										cout << "Error, while the building of " << MATRIX[0].tag << endl;
										state = 1;
										break;// To evit the incrementation of next matrix id
									};
								} else {
									temp2 = MATRIX[0].get_decimal(data);
									if (temp2 == 0 && data != "0") {
										state = 1;
										break;
									} else {
										MATRIX[0].build_from_scalar(temp2);
									}; 
								};
								// We save the operation
								MATRIX[MLS] = MATRIX[0];// the tag of the temp matrix will be restored in the end
								// I verify if it is a new matrix or not
								if (MLS==next_matrix_id) next_matrix_id++;
							} else {
								cout << "Error, the maximun number of matrices is attempt;" << endl;
								state = 1;
								break;
							};
							break;
						case 20:
							if ((temp=getMatrixIdByTag(MATRIX, data, next_matrix_id), temp) == MAX_NB_MATRIX && op != '.') {
								ERROR_MATRIX_NOT_FOUND(data);
								state = 1;
								break;
							} else if (op == '.' && (op=is_special_operator(data), op) == '\0') {
								cout << "Error: attribute " << data << " unknowed;";
								state = 1;
								break;
							} else {
								switch (op) {
									case '*':
										MATRIX[0] = MATRIX[MLS] * MATRIX[temp];
										break;
									case '+':
										MATRIX[0] = MATRIX[MLS] + MATRIX[temp];
										break;
									case '-':
										MATRIX[0] = MATRIX[MLS] - MATRIX[temp];
										break;
									case '/':
										MATRIX[0] = MATRIX[MLS] / MATRIX[temp];
										break;
									case '&':
										cout << (MATRIX[MLS] == MATRIX[temp]) << endl;
										break;
									case 'x':
										cout << MATRIX[MLS].trace() << endl;
										break;
									case 'd':
										cout << MATRIX[MLS].determinant() << endl;
										break;
									case 'a':
										MATRIX[0] = MATRIX[MLS].adjacent();
										break;
									case 'c':
										MATRIX[0] = MATRIX[MLS].cofactor();
										break;
									case 'i':
										MATRIX[0] = MATRIX[MLS].invert();
										break;
									case 't':
										MATRIX[0] = MATRIX[MLS].transpose();
										break;
									default:
										cout << "Internal error, operator;" << endl;
										state = 1;
								};
								if (state != 1) {
									if (op != '&' && op != 'x' && op != 'd') {
										// I print the result
										MATRIX[0].print();
									};
								};
							};
							break;
						default:
							cout << "Syntax error;" << endl;
					};
					// I clean data
					data = "";state = 0;
					// I continue getting instructions
					continue;
				} else if (cmd[current_char] == ':' && state == 0) {
					temp = getMatrixIdByTag(MATRIX, data, next_matrix_id);
					MLS = (temp != MAX_NB_MATRIX) ? temp : next_matrix_id;
					MATRIX[0].tag = data;
					data = "";
					state = 10;
				} else if (cmd[current_char] == '=' && state == 10) {
					// The order of temp matrix, don't change if not specified
					MATRIX[0].order = (data == "") ? MATRIX[0].order : short(MATRIX[0].get_decimal(data));
					if (MATRIX[0].order > MAX_ORDER) {
						cout << "Error: The maximun order is " << MAX_ORDER << ";" << endl;
						state = 1;
					} else {
						data = "";
						state = 11;
					};
				} else if (cmd[current_char] == '{' && state == 11) {
					data = "";
					state = 12;
				} else if (is_operator(cmd[current_char]) && state == 0) {
					op = cmd[current_char];
					MLS = getMatrixIdByTag(MATRIX, data, next_matrix_id);
					if (MLS != MAX_NB_MATRIX) {
						data = "";
						state = 20;
					} else {
						ERROR_MATRIX_NOT_FOUND(data);
						state = 1;
					};
				} else if (cmd[current_char] == '}' && state == 12) {
					state = 13;
				} else if (state == 12  || state == 11 || (is_character(cmd[current_char]) && state != 10) || (is_digit(cmd[current_char]) && state == 10)) {
					data += cmd[current_char];
				} else {
					cout << "Error, syntax error: <" << cmd[current_char] << "> unexpected" << endl;
					state = 1;
				};
				// I reassign the tag, because can be lost after some operations
				if (state == 0 || state == 1) {
					MATRIX[0].tag = '$';
				};
			};
			if (cmd == data && data != "") {
				// I set additionnal cmd
				if (cmd == "HELP") {
					cout << "No disponible for the moment" << endl;
				} else {
					temp = getMatrixIdByTag(MATRIX, data, next_matrix_id);
					if (temp != MAX_NB_MATRIX) {
						MATRIX[temp].print();
					} else {
						cout << "Error, keyword unknowed;" << endl;
						state = 1;
					};
				};
			};
			if (state == 1) {
				cout << "Line " << 1 << " and cell " << current_char << endl;
			};
			cout << ':' << state << endl;
		} else break;
	};
	cout << "GoodBye!" << endl;
    return 0;
};
