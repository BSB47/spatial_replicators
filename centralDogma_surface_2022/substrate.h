#ifndef SUBSTRATE_H
#define SUBSTRATE_H
class Substrate{
private:
	enum typeInComplex
	{
		free,
		occupied,
	};

	typeInComplex m_typeComp{};

	int m_xcor{};
	int m_ycor{};

	double m_mutation_probability{};
	
public: 
	Substrate(typeInComplex type=free, int xcor, int ycor) : m_typeComp{type}, m_xcor{xcor}, m_ycor{ycor}
	{
	}

	int getXcor(){return m_xcor;}
	int getYcor(){return m_ycor;}
};
#endif
