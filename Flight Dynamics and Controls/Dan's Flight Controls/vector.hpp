// Vector class - a few linear algebra functions for use in autopilot.

/**********************************************************
*    Structs for passing commands
*    - These are a couple of structs I used in writing my code; they're just proposed 
*    - versions until we can nail something down more firmly.
*
**********************************************************/

struct pathParams{  // Output parameters
	double VaC;
	double heightCommanded; // Positive or negative? 
	double chiCommandedRadians;
	double phiFFradians;
};

struct threeDVector{
	double xComponent;
	double yComponent;
	double zComponent;
};

struct pathInputParams{
	// Path planning information
	int orbitFlag; // One if by land, two if by sea... er, STRAIGHT for straight, ORBIT for orbit. 
	double VaD; // Desired Va
	threeDVector rPath;
	threeDVector qPath;
	threeDVector cOrbit;
	double orbitRadius;
	int orbitDirection; // CLOCKWISE or CNTCLOCKWISE
	// Current state information
	double pn;
	double pe;
	double height; // Height defined as negative? 
	double VaCurrent;
	// double alpha;
	// double beta;
	double phiRadians;
	double thetaRadians;
	double chiRadians;
	// double p;
	// double q;
	double r;
	// double Vg;
	// double wn
	// double we
	// double psi
	// Time? 
};

/**********************************************************
*    Class definitions
*    - Two class definitions; class Vector defines a 3x1
*    - vector, while class squareMatrix defines a 3x3 matrix.
*
**********************************************************/

class Vector{
private:
	// Arrays of the values in the vector and the normalized vector
	double* values;
	double* normalized;
public:
	Vector(double a, double b, double c){
		values = new double[3];
		normalized = new double[3];
		values[0] = a;
		values[1] = b;
		values[2] = c;
		normalized[0] = a / norm();
		normalized[1] = b / norm();
		normalized[2] = c / norm();
	}

	double norm(){
		// Return the induced norm, i.e. the square root of the inner product.
		return std::sqrt(pow(values[0], 2) + pow(values[1], 2) + pow(values[2], 2) );
	}
	double* getVector(){
		// Return an array of the values
		return values;
	}
	Vector getNormalized(){
		// Return an array of the normalized values
		return Vector(normalized[0], normalized[1], normalized[2]);
	}
	double operator[](int index){
		// Index operator. Return value[index].
		return values[index];
	}
	~Vector(){
		// Class destructor.
		delete[] values;
		delete[] normalized;
	}
};

class squareMatrix{
private:
	// Array of doubles containing the matrix values.
	double* values;
public:
	squareMatrix(double a11, double a12, double a13, double a21, double a22, double a23, double a31, double a32, double a33){
		// Initialization: First 3 values are first row, second 3 values are second row, third 3 values are third row. 
		values = new double[9];
		values[0] = a11;
		values[1] = a12;
		values[2] = a13;
		values[3] = a21;
		values[4] = a22;
		values[5] = a23;
		values[6] = a31;
		values[7] = a32;
		values[8] = a33;
	};
	double* getMatrix(){
		// Return an array of the values
		return values;
	}
	double operator[](int index){
		// Index operator. Returns the index in the array. I'm currently trying to figure out a [row][column] syntax as well.
		return values[index];
	}
/*	double operator[][](int row, int column){
		return 1.0;
	}*/
	squareMatrix transpose(){
		// Returns a square matrix that is the transpose of the original.
		return squareMatrix(values[0], values[3], values[6], values[1], values[4], values[7], values[2], values[5], values[8]);
	}
	~squareMatrix(){
		// Class destructor. 
		delete[] values;
	}
};


Vector operator*(squareMatrix a, Vector b){
	// Multiplication operator, square matrix times a vector. (3x3 times 3x1)
	double a1 = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	double a2 = a[3] * b[0] + a[4] * b[1] + a[5] * b[2];
	double a3 = a[6] * b[0] + a[7] * b[1] + a[8] * b[2];
	return Vector(a1, a2, a3);
}

Vector operator*(double scalar, Vector a){
	// Scalar multiply
	return Vector(scalar * a[0], scalar * a[1], scalar * a[2]);
}

Vector operator*(Vector a, double scalar){
	// Scalar multiply, opposite order.
	return Vector(scalar * a[0], scalar * a[1], scalar * a[2]);
}

Vector operator+(Vector a, Vector b){
	// Simple vector addition
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

Vector operator-(Vector a, Vector b){
	// Simple vector subtraction
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

squareMatrix operator*(squareMatrix a, squareMatrix b){
	// Square matrix multiplication, using square * column multiplication to reduce *my* work.
	Vector col1 = Vector(b[0], b[3], b[6]);
	Vector col2 = Vector(b[1], b[4], b[7]);
	Vector col3 = Vector(b[2], b[5], b[8]);

	Vector result1 = a * col1;
	Vector result2 = a * col2;
	Vector result3 = a * col3;  // Just re-using my vector multiply, since it still works.

	return squareMatrix(result1[0], result2[0], result3[0], result1[1], result2[1], result3[1], result1[2], result2[2], result3[2]);
}

double norm(Vector in){
	// Another norm function, just because
	return std::sqrt(pow(in[0],2) + pow(in[1], 2) + pow(in[2],2));
}

Vector normalize(Vector in){
	// Normalize a vector and return the normalized vector. Again, sort of a repeat, but it seems more useful here
	return in.getNormalized();
}

Vector crossProduct(Vector a, Vector b){
	// Compute the cross product of two vectors 
	double v0 = (a[1] * b[2]) - (a[2] * b[1]);
	double v1 = (a[2] * b[0]) - (a[0] * b[2]);
	double v2 = (a[0] * b[1]) - (a[1] * b[0]);

	return Vector(v0, v1, v2);
}

double dotProduct(Vector a, Vector b){
	// Compute the dot product of two vectors
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}