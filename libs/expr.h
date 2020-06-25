/**
 * simple LL(1) expression parser
 *
 * @author Tobias Weber <tweber@ill.fr>
 * @date 28-mar-20
 * @license GPLv3, see 'LICENSE' file
 * @desc The present version was forked on 28-Mar-2020 from the privately developed "misc" project (https://github.com/t-weber/misc).
 *
 * References:
 *	- https://de.wikipedia.org/wiki/LL(k)-Grammatik
 */

#ifndef __TLIBS2_EXPR_PARSER_H__
#define __TLIBS2_EXPR_PARSER_H__


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <regex>
#include <unordered_map>
#include <vector>
#include <memory>
#include <cmath>
#include <boost/math/special_functions/erf.hpp>

#include "units.h"


namespace tl2 {


template<typename t_num=double>
class ExprParser
{
public:
	t_num parse(const std::string& str)
	{
		m_istr = std::make_shared<std::istringstream>(str);
		next_lookahead();
		t_num result = plus_term();

		// check if there would be are more tokens available?
		next_lookahead();
		bool at_eof = (m_lookahead == (int)Token::TOK_INVALID || m_lookahead == (int)Token::TOK_END);
		if(!at_eof)
			throw std::underflow_error("Not all input tokens have been consumed.");

		return result;
	}


protected:
	// ------------------------------------------------------------------------
	// tables / functions
	// ------------------------------------------------------------------------

	// real functions with zero parameters
	template<class t_val,
		typename std::enable_if<std::is_floating_point<t_val>::value>::type* =nullptr>
	t_val call_func0(const std::string& strName)
	{
		static const std::unordered_map<std::string, t_val(*)()> s_funcs =
		{
		};

		return s_funcs.at(strName)();
	}


	// real functions with one parameter
	template<class t_val,
		typename std::enable_if<std::is_floating_point<t_val>::value>::type* =nullptr>
	t_val call_func1(const std::string& strName, t_val t)
	{
		static const std::unordered_map<std::string, t_val(*)(t_val)> s_funcs =
		{
			{ "sin", std::sin }, { "cos", std::cos }, { "tan", std::tan },
			{ "asin", std::asin }, { "acos", std::acos }, { "atan", std::atan },
			{ "sinh", std::sinh }, { "cosh", std::cosh }, { "tanh", std::tanh },
			{ "asinh", std::asinh }, { "acosh", std::acosh }, { "atanh", std::atanh },

			{ "sqrt", std::sqrt }, { "cbrt", std::cbrt },
			{ "exp", std::exp },
			{ "log", std::log }, { "log2", std::log2 }, { "log10", std::log10 },

			{ "erf", std::erf }, { "erfc", std::erfc }, { "erf_inv", boost::math::erf_inv },

			{ "round", std::round }, { "ceil", std::ceil }, { "floor", std::floor },
			{ "abs", std::abs },
		};

		return s_funcs.at(strName)(t);
	}


	// real functions with two parameters
	template<class t_val,
		typename std::enable_if<std::is_floating_point<t_val>::value>::type* =nullptr>
	t_val call_func2(const std::string& strName, t_val t1, t_val t2)
	{
		static const std::unordered_map<std::string, t_val(*)(t_val, t_val)> s_funcs =
		{
			{ "pow", std::pow }, { "atan2", std::atan2 },
			{ "mod", std::fmod },
		};

		return s_funcs.at(strName)(t1, t2);
	}


	// real constants
	template<class t_val,
		typename std::enable_if<std::is_floating_point<t_val>::value>::type* =nullptr>
	t_val get_const(const std::string& strName)
	{
		static const std::unordered_map<std::string, t_val> s_consts =
		{
			{ "pi", __pi<t_val> },
			{ "hbar",  t_val(hbar<t_val>/meV<t_val>/sec<t_val>) },	// hbar in [meV s]
			{ "kB",  t_val(kB<t_val>/meV<t_val>*kelvin<t_val>) },	// kB in [meV / K]
		};

		return s_consts.at(strName);
	}



	// alternative: int functions with zero parameters
	template<class t_val,
		typename std::enable_if<std::is_integral<t_val>::value>::type* =nullptr>
	t_val call_func0(const std::string& strName)
	{
		static const std::unordered_map<std::string, t_val(*)()> s_funcs =
		{
		};

		return s_funcs.at(strName)();
	}


	// alternative: int functions with one parameter
	template<class t_val,
		typename std::enable_if<std::is_integral<t_val>::value>::type* =nullptr>
	t_val call_func1(const std::string& strName, t_val t)
	{
		static const std::unordered_map<std::string, t_val(*)(t_val)> s_funcs =
		{
			{ "abs", std::abs },
		};

		return s_funcs.at(strName)(t);
	}


	// alternative: int functions with two parameters
	template<class t_val,
		typename std::enable_if<std::is_integral<t_val>::value>::type* =nullptr>
	t_val call_func2(const std::string& strName, t_val t1, t_val t2)
	{
		static const std::unordered_map<std::string, std::function<t_val(t_val, t_val)>> s_funcs =
		{
			{ "pow", [](t_val t1, t_val t2) -> t_val { return t_val(std::pow(t1, t2)); } },
			{ "mod", [](t_val t1, t_val t2) -> t_val { return t1%t2; } },
		};

		return s_funcs.at(strName)(t1, t2);
	}


	// alternative: int constants
	template<class t_val,
		typename std::enable_if<std::is_integral<t_val>::value>::type* =nullptr>
	t_val get_const(const std::string& strName)
	{
		/*static const std::unordered_map<std::string, t_val> s_consts =
		{
			{ "pi", 3 },
		};

		return s_consts.at(strName);*/

		throw std::out_of_range("Undefined constant.");
	}
	// ------------------------------------------------------------------------



	// ------------------------------------------------------------------------
	// Lexer
	// ------------------------------------------------------------------------
	enum class Token : int
	{
		TOK_NUM		= 1000,
		TOK_IDENT	= 1001,
		TOK_END		= 1002,

		TOK_INVALID	= 10000,
	};


	/**
	 * find all matching tokens for input string
	 */
	std::vector<std::pair<int, t_num>> get_matching_tokens(const std::string& str)
	{
		std::vector<std::pair<int, t_num>> matches;

		if constexpr(std::is_floating_point_v<t_num>)
		{	// real
			std::regex regex{"[0-9]+(\\.[0-9]*)?|\\.[0-9]+([eE][-+]?[0-9]+)?"};
			std::smatch smatch;
			if(std::regex_match(str, smatch, regex))
			{
				t_num val{};
				std::istringstream{str} >> val;
				matches.push_back(std::make_pair((int)Token::TOK_NUM, val));
			}
		}
		else if constexpr(std::is_integral_v<t_num>)
		{	// real
			std::regex regex{"[0-9]+"};
			std::smatch smatch;
			if(std::regex_match(str, smatch, regex))
			{
				t_num val{};
				std::istringstream{str} >> val;
				matches.push_back(std::make_pair((int)Token::TOK_NUM, val));
			}
		}
		else
		{
			throw std::invalid_argument("Invalid number type.");
		}

		{	// ident
			std::regex regex{"[A-Za-z_][A-Za-z0-9_]*"};
			std::smatch smatch;
			if(std::regex_match(str, smatch, regex))
				matches.push_back(std::make_pair((int)Token::TOK_IDENT, 0.));
		}

		{	// tokens represented by themselves
			if(str == "+" || str == "-" || str == "*" || str == "/" ||
				str == "%" || str == "^" || str == "(" || str == ")" || str == ",")
				matches.push_back(std::make_pair((int)str[0], 0.));
		}

		return matches;
	}


	/**
	 * @return [token, yylval, yytext]
	 */
	std::tuple<int, t_num, std::string> lex()
	{
		std::string input, longest_input;
		std::vector<std::pair<int, t_num>> longest_matching;

		// find longest matching token
		while(1)
		{
			char c = m_istr->get();

			if(m_istr->eof())
				break;
			// if outside any other match...
			if(longest_matching.size() == 0)
			{
				// ...ignore white spaces
				if(c==' ' || c=='\t')
					continue;
				// ...end on new line
				if(c=='\n')
					return std::make_tuple((int)Token::TOK_END, t_num{0}, longest_input);
			}

			input += c;
			auto matching = get_matching_tokens(input);
			if(matching.size())
			{
				longest_input = input;
				longest_matching = matching;

				if(m_istr->peek() == std::char_traits<char>::eof())
					break;
			}
			else
			{
				// no more matches
				m_istr->putback(c);
				break;
			}
		}

		// at EOF
		if(longest_matching.size() == 0 && (input.length() == 0 || m_istr->eof()))
		{
			return std::make_tuple((int)Token::TOK_END, t_num{0}, longest_input);
		}

		// nothing matches
		if(longest_matching.size() == 0)
		{
			std::ostringstream ostr;
			ostr << "Invalid input in lexer: \"" << input << "\".";
			throw std::runtime_error(ostr.str());
			return std::make_tuple((int)Token::TOK_INVALID, t_num{0}, longest_input);
		}

		// several possible matches
		if(longest_matching.size() > 1)
		{
			std::ostringstream ostr;
			ostr << "Warning: Ambiguous match in lexer for token \"" << longest_input << "\".";
			throw std::runtime_error(ostr.str());
		}

		// found match
		return std::make_tuple((int)std::get<0>(longest_matching[0]), std::get<1>(longest_matching[0]), longest_input);
	}
	// ------------------------------------------------------------------------



	// ----------------------------------------------------------------------------
	// Lexer interface
	// ----------------------------------------------------------------------------
	void next_lookahead()
	{
		std::tie(m_lookahead, m_lookahead_val, m_lookahead_text) = lex();
	}


	void match(int expected)
	{
		if(m_lookahead != expected)
		{
			std::ostringstream ostr;
			ostr << "Could not match symbol! Expected: " << expected << ", got: " << m_lookahead << ".";
			throw std::runtime_error(ostr.str());
		}
	}
	// ----------------------------------------------------------------------------



	// ----------------------------------------------------------------------------
	// Productions
	// ----------------------------------------------------------------------------
	/**
	 * +,- terms
	 * (lowest precedence, 1)
	 */
	t_num plus_term()
	{
		// plus_term -> mul_term plus_term_rest
		if(m_lookahead == '(' || m_lookahead == (int)Token::TOK_NUM || m_lookahead == (int)Token::TOK_IDENT)
		{
			t_num term_val = mul_term();
			t_num expr_rest_val = plus_term_rest(term_val);

			return expr_rest_val;
		}
		else if(m_lookahead == '+')	// unary +
		{
			next_lookahead();
			t_num term_val = mul_term();
			t_num expr_rest_val = plus_term_rest(term_val);

			return expr_rest_val;
		}
		else if(m_lookahead == '-')	// unary -
		{
			next_lookahead();
			t_num term_val = -mul_term();
			t_num expr_rest_val = plus_term_rest(term_val);

			return expr_rest_val;
		}

		if(m_lookahead == 0 || m_lookahead == EOF)
			exit(0);

		std::ostringstream ostr;
		ostr << "Invalid lookahead in " << __func__ << ": " << m_lookahead << ".";
		throw std::runtime_error(ostr.str());
		return 0.;
	}


	t_num plus_term_rest(t_num arg)
	{
		// plus_term_rest -> '+' mul_term plus_term_rest
		if(m_lookahead == '+')
		{
			next_lookahead();
			t_num term_val = arg + mul_term();
			t_num expr_rest_val = plus_term_rest(term_val);

			return expr_rest_val;
		}

		// plus_term_rest -> '-' mul_term plus_term_rest
		else if(m_lookahead == '-')
		{
			next_lookahead();
			t_num term_val = arg - mul_term();
			t_num expr_rest_val = plus_term_rest(term_val);

			return expr_rest_val;
		}
		// plus_term_rest -> epsilon
		else if(m_lookahead == ')' || m_lookahead == (int)Token::TOK_END || m_lookahead == ',')
		{
			return arg;
		}

		std::ostringstream ostr;
		ostr << "Invalid lookahead in " << __func__ << ": " << m_lookahead << ".";
		throw std::runtime_error(ostr.str());
		return 0.;
	}


	/**
	 * *,/,% terms
	 * (precedence 2)
	 */
	t_num mul_term()
	{
		// mul_term -> pow_term mul_term_rest
		if(m_lookahead == '(' || m_lookahead == (int)Token::TOK_NUM || m_lookahead == (int)Token::TOK_IDENT)
		{
			t_num factor_val = pow_term();
			t_num term_rest_val = mul_term_rest(factor_val);

			return term_rest_val;
		}

		std::ostringstream ostr;
		ostr << "Invalid lookahead in " << __func__ << ": " << m_lookahead << ".";
		throw std::runtime_error(ostr.str());
		return 0.;
	}


	t_num mul_term_rest(t_num arg)
	{
		// mul_term_rest -> '*' pow_term mul_term_rest
		if(m_lookahead == '*')
		{
			next_lookahead();
			t_num factor_val = arg * pow_term();
			t_num term_rest_val = mul_term_rest(factor_val);

			return term_rest_val;
		}

		// mul_term_rest -> '/' pow_term mul_term_rest
		else if(m_lookahead == '/')
		{
			next_lookahead();
			t_num factor_val = arg / pow_term();
			t_num term_rest_val = mul_term_rest(factor_val);

			return term_rest_val;
		}

		// mul_term_rest -> '%' pow_term mul_term_rest
		else if(m_lookahead == '%')
		{
			next_lookahead();
			t_num factor_val = std::fmod(arg, pow_term());
			t_num term_rest_val = mul_term_rest(factor_val);

			return term_rest_val;
		}

		// mul_term_rest -> epsilon
		else if(m_lookahead == '+' || m_lookahead == '-' || m_lookahead == ')'
			|| m_lookahead == (int)Token::TOK_END || m_lookahead == ',')
		{
			return arg;
		}

		std::ostringstream ostr;
		ostr << "Invalid lookahead in " << __func__ << ": " << m_lookahead << ".";
		throw std::runtime_error(ostr.str());
		return 0.;
	}


	/**
	 * ^ terms
	 * (precedence 3)
	 */
	t_num pow_term()
	{
		// pow_term -> factor pow_term_rest
		if(m_lookahead == '(' || m_lookahead == (int)Token::TOK_NUM || m_lookahead == (int)Token::TOK_IDENT)
		{
			t_num factor_val = factor();
			t_num term_rest_val = pow_term_rest(factor_val);

			return term_rest_val;
		}

		std::ostringstream ostr;
		ostr << "Invalid lookahead in " << __func__ << ": " << m_lookahead << ".";
		throw std::runtime_error(ostr.str());
		return 0.;
	}


	t_num pow_term_rest(t_num arg)
	{
		// pow_term_rest -> '^' factor pow_term_rest
		if(m_lookahead == '^')
		{
			next_lookahead();
			t_num factor_val = std::pow(arg, factor());
			t_num term_rest_val = pow_term_rest(factor_val);

			return term_rest_val;
		}

		// pow_term_rest -> epsilon
		else if(m_lookahead == '+' || m_lookahead == '-' || m_lookahead == ')'
			|| m_lookahead == (int)Token::TOK_END || m_lookahead == ','
			|| m_lookahead == '*' || m_lookahead == '/' || m_lookahead == '%')
		{
			return arg;
		}

		std::ostringstream ostr;
		ostr << "Invalid lookahead in " << __func__ << ": " << m_lookahead << ".";
		throw std::runtime_error(ostr.str());
		return 0.;
	}


	/**
	 * () terms, real factor or identifier
	 * (highest precedence, 4)
	 */
	t_num factor()
	{
		// factor -> '(' plus_term ')'
		if(m_lookahead == '(')
		{
			next_lookahead();
			t_num expr_val = plus_term();
			match(')');
			next_lookahead();

			return expr_val;
		}

		// factor -> TOK_NUM
		else if(m_lookahead == (int)Token::TOK_NUM)
		{
			t_num val = m_lookahead_val;
			next_lookahead();

			return val;
		}

		// factor -> TOK_IDENT
		else if(m_lookahead == (int)Token::TOK_IDENT)
		{
			const std::string ident = m_lookahead_text;
			next_lookahead();

			// function call
			// using next m_lookahead, grammar still ll(1)?
			if(m_lookahead == '(')
			{
				next_lookahead();

				// 0-argument function
				// factor -> TOK_IDENT '(' ')'
				if(m_lookahead == ')')
				{
					next_lookahead();

					return call_func0<t_num>(ident);
				}

				// function with arguments
				else
				{
					// first argument
					t_num expr_val1 = plus_term();

					// one-argument-function
					// factor -> TOK_IDENT '(' plus_term ')'
					if(m_lookahead == ')')
					{
						next_lookahead();

						return call_func1<t_num>(ident, expr_val1);
					}

					// two-argument-function
					// factor -> TOK_IDENT '(' plus_term ',' plus_term ')'
					else if(m_lookahead == ',')
					{
						next_lookahead();
						t_num expr_val2 = plus_term();
						match(')');
						next_lookahead();

						return call_func2<t_num>(ident, expr_val1, expr_val2);
					}
					else
					{
						std::ostringstream ostr;
						ostr << "Invalid function call to \"" << ident << "\".";
						throw std::runtime_error(ostr.str());
					}
				}
			}

			// variable lookup
			else
			{
				return get_const<t_num>(ident);
			}
		}

		std::ostringstream ostr;
		ostr << "Invalid lookahead in " << __func__ << ": " << m_lookahead << ".";
		throw std::runtime_error(ostr.str());
		return 0.;
	}
	// ----------------------------------------------------------------------------


private:
	std::shared_ptr<std::istream> m_istr;

	int m_lookahead = (int)Token::TOK_INVALID;
	t_num m_lookahead_val = 0;
	std::string m_lookahead_text;
};


}
#endif
