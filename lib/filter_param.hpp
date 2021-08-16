
/*
 * filter_param.hpp
 *
 *  Created on: 2021/07/22
 *      Author: matsu
 *
 * This cord is written by UTF-8
 */

#ifndef FILTER_PARAM_HPP_
#define FILTER_PARAM_HPP_

#define _USE_MATH_DEFINES

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <regex>
#include <complex>
#include <functional>
#include <random>

using namespace std;

template <typename... Args>
string format(const string &, Args...);

/* フィルタの種類を示す列挙体
 *
 */
enum class FilterType
{
	LPF,
	HPF,
	BPF,
	BEF,
	Other
};

/* バンド(周波数帯域)の種別を示す列挙体
 *   Pass : 通過域
 *   Stop : 阻止域
 *   Transition : 遷移域
 */
enum class BandType
{
	Pass,
	Stop,
	Transition
};

/* バンド(周波数帯域)の情報をまとめた構造体
 *   type : 帯域の種類(通過・阻止・遷移)
 *   left : 帯域の左端正規化周波数 [0:0.5)
 *   right : 帯域の右端正規化周波数 (0:0.5]
 *
 *   left < rightであること。それ以外の場合、エラー終了。
 */
struct BandParam
{
protected:
	BandType band_type;
	double left_side;
	double right_side;

public:
	BandParam(BandType input_type, double input_left, double input_right)
		: band_type(input_type), left_side(0.0), right_side(0.0)
	{
		// バンドの条件に合わないとき、エラー終了
		if (input_left < 0.0 || input_left > input_right || input_right > 0.5)
		{
			fprintf(stderr,
				"Error: [%s l.%d]Band edge is illegal(left :%6.3f, right :%6.3f)\n",
				__FILE__, __LINE__, input_left, input_right);
			exit(EXIT_FAILURE);
		}

		left_side = input_left;
		right_side = input_right;
	}

	BandType type() const
	{ return band_type; }
	double left() const
	{ return left_side; }
	double right() const
	{ return right_side; }
	double width() const
	{ return right_side - left_side; }
	string sprint();
};

struct FilterParam
{
protected:
	// フィルタパラメータ

	unsigned int n_order;
	unsigned int m_order;
	vector<BandParam> bands;
	unsigned int nsplit_approx;
	unsigned int nsplit_transition;
	double group_delay;
	double threshold_riple;

	// 内部パラメータ
	
	vector<vector<complex<double>>> csw;			// 複素正弦波e^-jωを周波数帯域別に格納
	vector<vector<complex<double>>> csw2;			// 複素正弦波e^-j2ωを周波数帯域別に格納
	vector<vector<complex<double>>> desire_res;		// 所望特性の周波数特性

	function< vector<vector<complex<double>>>(const FilterParam*, const vector<double>&) > freq_res_func;
	function< double(const FilterParam*, const vector<double>&) > stability_func;

	// 内部メソッド

	FilterParam()
	:n_order(0), m_order(0),
	 nsplit_approx(0), nsplit_transition(0), group_delay(0.0),
	 threshold_riple(1.0)
	{}
	vector<vector<complex<double>>> freq_res_se(const vector<double>&) const;
	vector<vector<complex<double>>> freq_res_no(const vector<double>&) const;
	vector<vector<complex<double>>> freq_res_mo(const vector<double>&) const;


	double judge_stability_even(const vector<double>&) const;
	double judge_stability_odd(const vector<double>&) const;

public:
	FilterParam(unsigned int, unsigned int, vector<BandParam>,
				unsigned int, unsigned int, double);



	// get function

	unsigned int pole_order() const
	{ return m_order; }
	unsigned int zero_order() const
	{ return n_order; }
	unsigned int opt_order() const
	{ return 1 + n_order + m_order; }
	vector<BandParam> fbands() const
	{ return bands; }
	unsigned int partition_approx() const
	{ return nsplit_approx; }
	unsigned int partition_transition() const
	{ return nsplit_transition; }
	double gd() const
	{ return group_delay; }

	// set function
	/* # フィルタ構造体
	 *   振幅隆起として検知する閾値を変更する
	 *   デフォルト値は1.0
	 */
	void set_threshold_riple(double input)
	{ threshold_riple = input; }

	// normal function
	/* フィルタ構造体
	 *   周波数特性計算関数
	 *   コンストラクタに与えられた周波数帯域に
	 *   応じて，縦続型IIRフィルタの周波数特性を計算する
	 *   また，係数列も次数によって適宜分割される
	 *
	 *   # 引数
	 *   vector<double> coef : 係数列(a0, a1, a2[0], a2[1],..., b1, b2[0], b2[1],...)
	 *   #返り値
	 *   vector<vector<complex<double>>> response : 周波数帯域-周波数分割数の2重配列
	 */
	vector<vector<complex<double>>> freq_res(const vector<double>& coef) const
	{ return this->freq_res_func(this, coef); }
	
	/*
	 *
	 */
	double judge_stability(const vector<double>& coef) const
	{ return this->stability_func(this, coef); }

	double evaluate(const vector<double>&) const;
	vector<double> init_coef(const double, const double, const double) const;
	vector<double> init_stable_coef(const double, const double) const;

	// static function
	
	static vector<FilterParam> read_csv(string&);

	template <typename... Args>
	static vector<BandParam> gen_bands(FilterType, Args...);
	static FilterType analyze_type(const string&);
	static vector<double> analyze_edges(const string&);
	static vector<complex<double>> gen_csw(const BandParam&, const unsigned int);
	static vector<complex<double>> gen_csw2(const BandParam&, const unsigned int);
	static vector<complex<double>> gen_desire_res(const BandParam&, const unsigned int, const double);
};

//-------template function---------------------------------------
/* # String format function
 *
 */
template <typename... Args>
string format(const string &fmt, Args... args)
{
	size_t len = snprintf(nullptr, 0, fmt.c_str(), args...);
	vector<char> buf(len + 1);
	snprintf(&buf[0], len + 1, fmt.c_str(), args...);
	return string(&buf[0], &buf[0] + len);
}

/* # フィルタ構造体
 *   フィルタのタイプと周波数帯域端により，
 *   連続した周波数帯域の配列を生成する
 *
 * # 引数
 * FilterType ftype : フィルタタイプ(L.P.F.やH.P.F.など)
 * Args... edges : doubleの可変長引数を想定
 *                   フィルタタイプに合わない帯域端，
 *                   フィルタタイプが`Other`の場合には，
 *                   異常終了する
 * # 返り値
 * vector<BandParam> bands : 連続した周波数帯域の配列
 */
template <typename... Args>
vector<BandParam> FilterParam::gen_bands(FilterType ftype, Args... edges)
{
	vector<BandParam> bands;

	switch (ftype)
	{
		case FilterType::LPF:
		{
			int nedge = sizeof...(edges);
			if (nedge != 2)
			{
				fprintf(stderr,
					"Error: [%s l.%d]It has not been implement yet.\n",
					__FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}

			double edge[] = {edges...};
			bands.emplace_back(BandParam(BandType::Pass, 0.0, edge[0]));
			bands.emplace_back(BandParam(BandType::Transition, edge[0], edge[1]));
			bands.emplace_back(BandParam(BandType::Stop, edge[1], 0.5));
			break;
		}
		case FilterType::Other:
		{
			fprintf(stderr, 
				"Error: [%s l.%d]Other type filter can't use this function. "
				"Please, make `vector<BandParam>` your self.\n",
				__FILE__, __LINE__);
			exit(EXIT_FAILURE);
			break;
		}
		default:
		{
			fprintf(stderr,
				"Error: [%s l.%d]It has not been implement yet.\n",
				__FILE__, __LINE__);
			exit(EXIT_FAILURE);
			break;
		}
	}

	return bands;
}

#endif /* FILTER_PARAM_HPP_ */
