/*
 * filter_param.cpp
 *
 *  Created on: 2021/07/22
 *      Author: matsu
 */

#include "filter_param.hpp"

using namespace std;

FILE *fileopen(const string &filename, const char mode, const string &call_file, const int call_line)
{
	FILE *fp = fopen(filename.c_str(), &mode);
	if (fp == NULL)
	{
		fprintf(stderr,
				"Error: [%s l.%d]Can't open file.(file name : %s, mode : %c)\n",
				call_file.c_str(), call_line, filename.c_str(), mode);
		exit(EXIT_FAILURE);
	}
	return fp;
}

vector<string> split_char(string &str, char spliter)
{
	stringstream ss{str};
	vector<std::string> strs;
	string buf;
	while (std::getline(ss, buf, spliter))
	{
		strs.emplace_back(buf);
	}
	return strs;
}

string BandParam::sprint()
{
	string str;
	switch (band_type)
	{
	case BandType::Pass:
	{
		auto fmt = string("PassBand(%6.3f, %6.3f)");
		str = format(fmt, left_side, right_side);
		break;
	}
	case BandType::Stop:
	{
		str = format("StopBand(%6.3f, %6.3f)", left_side, right_side);
		break;
	}
	case BandType::Transition:
	{
		str = format("TransitionBand(%6.3f, %6.3f)", left_side, right_side);
		break;
	}
	}
	return str;
}

FilterParam::FilterParam
(unsigned int zero, unsigned int pole, vector<BandParam> input_bands,
		unsigned int input_nsplit_approx, unsigned int input_nsplit_transition,
		double gd)
:n_order(zero), m_order(pole), bands(input_bands),
 nsplit_approx(input_nsplit_approx), nsplit_transition(input_nsplit_transition),
 group_delay(gd),
 threshold_riple(1.0)
{
	// 周波数帯域の整合性チェック
	double band_left = 0.0;
	for (auto bp : bands)
	{
		if (bp.left() != band_left)
		{
			fprintf(stderr, "Error: [%s l.%d]Adjacent band edge is necessary same.\n"
							"Or first band of left side must be 0.0.\n",
					__FILE__, __LINE__);
			for (auto bp : bands)
			{
				fprintf(stderr, "%s\n", bp.sprint().c_str());
			}
			exit(EXIT_FAILURE);
		}
		band_left = bp.right();
	}
	if (bands.back().right() != 0.5)
	{
		fprintf(stderr,
				"Error: [%s l.%d]Last band of right side must be 0.5.\n",
				__FILE__, __LINE__);
		for (auto bp : bands)
		{
			fprintf(stderr, "%s\n", bp.sprint().c_str());
		}
		exit(EXIT_FAILURE);
	}

	// 帯域ごとの分割数算出
	double approx_range = 0.0;
	double transition_range = 0.0;

	for (auto bp : bands)
	{
		switch (bp.type())
		{
		case BandType::Pass:
		{
			approx_range += bp.width();
			break;
		}
		case BandType::Stop:
		{
			approx_range += bp.width();
			break;
		}
		case BandType::Transition:
		{
			transition_range += bp.width();
			break;
		}
		}
	}

	vector<unsigned int> split;
	split.reserve(bands.size());
	for (auto bp : bands)
	{
		switch (bp.type())
		{
		case BandType::Pass:
		{
			split.emplace_back(
				(unsigned int)((double)nsplit_approx * bp.width() / approx_range));
			break;
		}
		case BandType::Stop:
		{
			split.emplace_back(
				(unsigned int)((double)nsplit_approx * bp.width() / approx_range));
			break;
		}
		case BandType::Transition:
		{
			split.emplace_back(
				(unsigned int)((double)nsplit_transition * bp.width() / transition_range));
			break;
		}
		}
	}
	split.at(0) += 1;

	// generate complex sin wave(e^-jω)
	csw.reserve(bands.size());
	csw2.reserve(bands.size());
	for (unsigned int i = 0; i < bands.size(); ++i)
	{
		csw.emplace_back(FilterParam::gen_csw(bands.at(i), split.at(i)));
		csw2.emplace_back(FilterParam::gen_csw2(bands.at(i), split.at(i)));
	}

	// decide using function
	if ((n_order % 2) == 0)
	{
		if ((m_order % 2) == 0)
		{
			freq_res_func = bind(&FilterParam::freq_res_se, this,placeholders::_1);
		}
		else
		{
			freq_res_func = bind(&FilterParam::freq_res_mo, this,placeholders::_1);
		}
	}
	else
	{
		if ((m_order % 2) == 0)
		{
			freq_res_func = bind(&FilterParam::freq_res_no, this,placeholders::_1);
		}
		else
		{
			printf("so is unimplemented.\n");
		}
	}
}

vector<FilterParam> FilterParam::read_csv(string &filename)
{
	vector<FilterParam> filter_params;
	ifstream ifs(filename);
	if (!ifs)
	{
		fprintf(stderr,
				"Error: [%s l.%d]Can't open file.(file name : %s, mode : read)\n",
				__FILE__, __LINE__, filename.c_str());
		exit(EXIT_FAILURE);
	}

	string buf;
	getline(ifs, buf); // ヘッダ読み飛ばし
	while (getline(ifs, buf))
	{
		auto vals = split_char(buf, ',');
		auto type = FilterParam::analyze_type(vals.at(3));
		vector<BandParam> bands;

		switch (type)
		{
		case FilterType::LPF:
		{
			auto edges = FilterParam::analyze_edges(vals.at(3));
			if (edges.size() != 2)
			{
				fprintf(stderr, "Error: [%s l.%d]Format of filter state is illegal.(input : \"%s\")\n"
								"If you assign filter type L.P.F. , length of edges is only 2.\n",
						__FILE__, __LINE__, vals.at(3).c_str());
				exit(EXIT_FAILURE);
			}
			bands = FilterParam::gen_bands(FilterType::LPF, edges.at(0), edges.at(1));
			break;
		}
		default:
		{
			fprintf(stderr, "Error: [%s l.%d]It has not been implement yet.(input : \"%s\")\n",
					__FILE__, __LINE__, vals.at(3).c_str());
			exit(EXIT_FAILURE);
		}
		}

		filter_params.emplace_back(
				FilterParam(atoi(vals.at(1).c_str()), atoi(vals.at(2).c_str()),
						bands, atoi(vals.at(5).c_str()),
						atoi(vals.at(6).c_str()), atof(vals.at(4).c_str())));
	}

	// 次数に応じた関数を格納

	return filter_params;
}

/* フィルタのタイプを簡易入力する
 *   LPF(0.2, 0.3) : L.P.F.で通過域端 0.2, 阻止域端 0.3のフィルタ
 *
 *   Other(path) : その他の特性のフィルタ。pathは周波数帯域の情報を記したファイルへの相対パス
 */
FilterType FilterParam::analyze_type(string &input)
{
	FilterType type;
	string str = input;
	str.erase(remove(str.begin(), str.end(), ' '), str.end());
	string csv_space = regex_replace(str, regex("[(),:]+"), " ");
	auto input_type = split_char(csv_space, ' ');

	if (input_type.size() <= 1)
	{
		fprintf(stderr,
				"Error: [%s l.%d]Format of filter type is illegal.(input : \"%s\")\n",
				__FILE__, __LINE__, input.c_str());
		exit(EXIT_FAILURE);
	}
	else if (input_type.at(0) == "LPF")
	{
		type = FilterType::LPF;
	}
	else if (true)	// hpf, bpf, bef バリエーション
	{
		fprintf(stderr,
				"Error: [%s l.%d]It has not been implement yet.(input : \"%s\")\n",
				__FILE__, __LINE__, input.c_str());
		exit(EXIT_FAILURE);
	}
	else
	{
		fprintf(stderr,
				"Error: [%s l.%d]Filter Type is undefined.(input : \"%s\")\n",
				__FILE__, __LINE__, input.c_str());
		exit(EXIT_FAILURE);
	}

	return type;
}

vector<double> FilterParam::analyze_edges(string &input)
{
	vector<double> edge;
	string str = input;
	str.erase(remove(str.begin(), str.end(), ' '), str.end());
	string csv_space = regex_replace(str, regex("[(),:]+"), " ");
	auto input_type = split_char(csv_space, ' ');
	input_type.erase(input_type.begin());

	if (input_type.size() <= 0)
	{
		fprintf(stderr,
				"Error: [%s l.%d]Format of filter state is illegal.(input : \"%s\")\n",
				__FILE__, __LINE__, input.c_str());
		exit(EXIT_FAILURE);
	}
	for (auto v : input_type)
	{
		edge.emplace_back(atof(v.c_str()));
	}
	return edge;
}

vector<complex<double>> FilterParam::gen_csw(BandParam &bp, unsigned int nsplit)
{
	vector<complex<double>> csw;
	csw.reserve(nsplit);
	double step_size = bp.width() / nsplit;
	double left = bp.left();
	const double dpi = -2.0 * M_PI;			// double pi

	for (unsigned int i = 0; i < nsplit; ++i)
	{
		csw.emplace_back(polar(1.0, dpi * (left + step_size * (double) i)));
	}

	return csw;
}

vector<complex<double>> FilterParam::gen_csw2(BandParam &bp,
		unsigned int nsplit)
{
	vector<complex<double>> csw2;
	csw2.reserve(nsplit);
	double step_size = bp.width() / nsplit;
	double left = bp.left();
	const double dpi = -4.0 * M_PI;			// quadrical pi

	for (unsigned int i = 0; i < nsplit; ++i)
	{
		csw2.emplace_back(polar(1.0, dpi * (left + step_size * (double) i)));
	}

	return csw2;
}

vector<vector<complex<double>>> FilterParam::freq_res_no(vector<double>& coef)    // 周波数特性計算関数
{
	vector<vector<complex<double>>> freq;
	vector<vector<complex<double>>> freq_denominator;
	vector<vector<complex<double>>> freq_numerator;
	complex<double> one(1.0 , 0.0);

    for(unsigned int i = 0; i < bands.size() ; ++i)    // 周波数帯域のループ(L.P.F.なら３つ)
    {
      // csw.at(i), csw2.at(i), bands.at(i)がその周波数帯域で使う値に
      // cswは複素正弦波、e^-jω
      
      for(unsigned int j = 0; j < csw.at(i).size() ; ++j)    // 周波数帯域内の分割数によるループ
      {
		//coef[0]=a_0,coef[1]=a_1,coef[2]=b_1,coef[3]=a_2,coef[4]=b_2...etc
		//奇数はa,偶数はb

	    freq_numerator.at(i).at(j) = one + coef[1] * csw.at(i).at(j);
		  for(unsigned int N = 2; N < coef.size() ; ++N)	//フィルタ係数の計算用ループ(分子)
		  {
			if(N % 2 == 0)
			freq_numerator.at(i).at(j) = freq_numerator.at(i).at(j) * (one + coef[2 * N - 1] * csw.at(i).at(j) + coef[2 * N + 1] * csw2.at(i).at(j));
		  }

		freq_denominator.at(i).at(j) = one + coef[2] * csw.at(i).at(j);
		  for(unsigned int M = 2; M < coef.size() ; ++M)	//フィルタ係数の計算用ループ(分母)
		  {
			if(M % 2 == 0)
			freq_denominator.at(i).at(j) = freq_denominator.at(i).at(j) * (one + coef[2 * M] * csw.at(i).at(j) + coef[2 * M + 2] * csw2.at(i).at(j));
		  }

			freq.at(i).at(j) = coef[0] * freq_numerator.at(i).at(j) / freq_denominator.at(i).at(j);

		}
	 } 
  return freq ;
}
vector<vector<complex<double>>> FilterParam::freq_res_se(vector<double>& coef)
{
	vector<vector<complex<double>>> res;
		res.reserve(bands.size());
	const complex<double> one(1.0, 0.0);

	for (unsigned int i = 0; i < bands.size(); ++i)    // 周波数帯域のループ
	{
		vector<complex<double>> band_res;
			band_res.reserve(csw.at(i).size());

		for (unsigned int j = 0; j < csw.at(i).size(); ++j)  // 周波数帯域内の分割数によるループ
		{
			complex<double> frac_over(1.0, 1.0);
			complex<double> frac_under(1.0, 1.0);

			for (unsigned int n = 1; n < n_order; n += 2)
			{
				frac_over *= one + coef.at(n)*csw.at(i).at(j) + coef.at(n + 1)*csw2.at(i).at(j);
			}
			for (unsigned int m = n_order + 1; m < opt_order(); m += 2)
			{
				frac_under *= one + coef.at(m)*csw.at(i).at(j) + coef.at(m + 1)*csw2.at(i).at(j);
			}
			band_res.emplace_back(coef.at(0) * (frac_over / frac_under));
		}
		res.emplace_back(band_res);
	}

	return res;
}

vector<vector<complex<double>>> FilterParam::freq_res_no(vector<double> &coef)
{
	vector<vector<complex<double>>> freq;
		freq.reserve(bands.size());
	const complex<double> one(1.0, 0.0);

	for (unsigned int i = 0; i < bands.size(); ++i)    // 周波数帯域のループ
	{
		vector<complex<double>> freq_band;
			freq_band.reserve(csw.at(i).size());

		for (unsigned int j = 0; j < csw.at(i).size(); ++j)  // 周波数帯域内の分割数によるループ
		{
			complex<double> freq_denominator(1.0, 1.0);
			complex<double> freq_numerator(1.0, 1.0);

			freq_numerator *= one + (coef.at(1) * csw.at(i).at(j));
			for (unsigned int n = 2; n < n_order; n += 2)		//分子の総乗ループ
			{
				freq_numerator *= (one + coef.at(n)*csw.at(i).at(j) + coef.at(n + 1)*csw2.at(i).at(j));
			}
			for (unsigned int m = n_order + 1; m < opt_order(); m += 2)		//分母の総乗ループ
			{
				freq_denominator *= (one + coef.at(m)*csw.at(i).at(j) + coef.at(m + 1)*csw2.at(i).at(j));
			}

			freq_band.emplace_back( coef.at(0) * (freq_numerator / freq_denominator));
		}
		freq.emplace_back(freq_band);
	}
	
	return freq;
}

vector<vector<complex<double>>> FilterParam::freq_res_mo(vector<double> &coef) // 周波数特性計算関数
{
	vector<vector<complex<double>>> freq;
	freq.reserve(bands.size());

	for (unsigned int i = 0; i < bands.size(); ++i) // 周波数帯域のループ(L.P.F.なら３つ)
	{
		// csw.at(i), csw2.at(i), bands.at(i)がその周波数帯域で使う値に
		// cswは複素正弦波、e^-jω

		vector<complex<double>> freq_band;
		freq_band.reserve(csw.at(i).size());

		for (unsigned int j = 0; j < csw.at(i).size(); ++j) // 周波数帯域内の分割数によるループ
		{
			// 2次の分子なら
			// 1 + coef[0]*csw.at(i).at(j) + coef[1]*csw2.at(i).at(j)
			// みたいにかける
			// 係数のインデックスはおかしいけど、適当に埋めてあるだけです
			complex<double> nume(1.0, 1.0);

			for (unsigned int n = 1; n < n_order; n += 2)
			{
				nume *= 1.0 + coef[n] * csw.at(i).at(j) + coef[n + 1] * csw2.at(i).at(j);
			}

			complex<double> deno(1.0, 1.0);

			deno *= 1.0 + (coef.at(n_order + 1) * csw.at(i).at(j));

			for (unsigned int m = n_order + 2; m < opt_order(); m += 2) //分母だけ奇数
			{
				deno *= 1.0 + coef[m] * csw.at(i).at(j) + coef[m + 1] * csw2.at(i).at(j);
			}

			freq_band.emplace_back(coef.at(0) * nume / deno);
		}
		freq.emplace_back(freq_band);
	}
	
	return freq;
}
