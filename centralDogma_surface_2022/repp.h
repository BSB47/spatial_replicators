#ifndef REPP_H
#define REPP_H
class Repp{
private:
	double m_kppp{}; //last two letters stand for template/product
	double m_kppq{};
	double m_kpqq{};
	double m_kpqp{};
	double m_kqpp{};
	double m_kqpq{};
	double m_kqqq{};
	double m_kqqp{};

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
	Repp() = default;

	Repp(double kppp=1, double kppq=1, double kpqq=1, double kpqp=1, double kqpp=1, double kqpq=1, double kqqq=1, double kqqp=1, typeInComplex type=free, int xcor, int ycor) : m_typeComp{type}, m_xcor{xcor}, m_ycor{ycor}
	{
	}

	int getXcor(){return m_xcor;}
	int getYcor(){return m_ycor;}
};
#endif
