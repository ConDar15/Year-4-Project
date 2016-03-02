
double geometric_acos(double x, unsigned int n)
{
	//Ensures the given value is a valid cosine value
	assert(x >= -1 && x <= 1);
	
	//Reversing the last line of geometric_cos_bounded
	double h = 2*(1-x);

	//Reverses the iterative process oof geometric_cos_bounded
	for(int i = 0; i < n; i++)
		h = 2 - sqrt(4 - h);

	//Reverses the initialisation proceduce in geometric_cos_bounded
	h *= pow(4, n);
	return sqrt(h);
}

double geometric_asin(double x, unsigned int n)
{
	assert(x >= -1 && x <= 1);
	return HALF_PI - geometric_acos(x, n);
}

double geometric_atan(double x, unsigned int n)
{
	return geometric_asin(x/sqrt(x*x + 1), n);
}
