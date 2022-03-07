#ifndef REPQ_H
#define REPQ_H
class Repq{
private:
	double m_kqpp{}; //last two letters stand for template/product
	double m_kqpq{};
	double m_kqqq{};
	double m_kqqp{};
	double m_kppp{}; 
	double m_kppq{};
	double m_kpqq{};
	double m_kpqp{};

	enum typeInComplex
	{
		free,
		temp,
		cata
	};

	typeInComplex m_typeComp{};

	int m_xcor{};
	int m_ycor{};

	double m_mutation_probability{};
	
public: 
	Repq() = default;

	Repq(double kppp=1, double kppq=1, double kpqq=1, double kpqp=1, double kqpp=1, double kqpq=1, double kqqq=1, double kqqp=1, typeInComplex type=free, int xcor, int ycor) : m_typeComp{type}, m_xcor{xcor}, m_ycor{ycor}
	{
	}

	int getXcor(){return m_xcor;}
	int getYcor(){return m_ycor;}
};
#endif
