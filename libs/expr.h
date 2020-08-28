/**
 * tlibs2 -- simple LL(1) expression parser
 * @author Tobias Weber <tweber@ill.fr>
 * @date 28-mar-2020
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

#if __has_include(<boost/math/special_functions/erf.hpp>)
	#include <boost/math/special_functions/erf.hpp>
	#define HAS_ERFINV
#endif

#include "units.h"


namespace tl2 {


template<class t_num=double>
t_num modfunc(t_num t1, t_num t2)
{
	if constexpr(std::is_floating_point_v<t_num>)
		return std::fmod(t1, t2);
	else if constexpr(std::is_integral_v<t_num>)
		return t1%t2;
	else
		throw std::runtime_error{"Invalid type for mod function."};
}



template<typename t_num=double> class ExprParser;


// ------------------------------------------------------------------------
// ast
// ------------------------------------------------------------------------
template<typename t_num=double>
class ExprAST
{
public:
	virtual ~ExprAST() = default;

	virtual t_num eval(const ExprParser<t_num>&) const = 0;
	virtual void print(std::ostream& ostr=std::cout, std::size_t indent=0) const = 0;
};


template<typename t_num=double>
class ExprASTBinOp : public ExprAST<t_num>
{
public:
	ExprASTBinOp(char op, const std::shared_ptr<ExprAST<t_num>>& left, const std::shared_ptr<ExprAST<t_num>>& right)
		: m_op{op}, m_left{left}, m_right{right}
	{}

	ExprASTBinOp(const ExprASTBinOp<t_num>&) = delete;
	const ExprASTBinOp<t_num>& operator=(const ExprASTBinOp<t_num>&) = delete;

	virtual ~ExprASTBinOp() = default;

	virtual t_num eval(const ExprParser<t_num>& context) const override
	{
		const t_num val_left = m_left->eval(context);
		const t_num val_right = m_right->eval(context);

		switch(m_op)
		{
			case '+': return val_left + val_right;
			case '-': return val_left - val_right;
			case '*': return val_left * val_right;
			case '/': return val_left / val_right;
			case '%': return modfunc<t_num>(val_left, val_right);
			case '^': return static_cast<t_num>(std::pow(val_left, val_right));
		}

		throw std::runtime_error{"Invalid binary operator."};
	}

	virtual void print(std::ostream& ostr=std::cout, std::size_t indent=0) const override
	{
		for(std::size_t i=0; i<indent; ++i) ostr << " | ";
		ostr << "binary operator " << m_op << "\n";
		m_left->print(ostr, indent+1);
		m_right->print(ostr, indent+1);
	}

private:
	char m_op{'+'};
	std::shared_ptr<ExprAST<t_num>> m_left{}, m_right{};
};


template<typename t_num=double>
class ExprASTUnOp : public ExprAST<t_num>
{
public:
	ExprASTUnOp(char op, const std::shared_ptr<ExprAST<t_num>>& child)
		: m_op{op}, m_child{child}
	{}

	ExprASTUnOp(const ExprASTUnOp<t_num>&) = delete;
	const ExprASTUnOp<t_num>& operator=(const ExprASTUnOp<t_num>&) = delete;

	virtual ~ExprASTUnOp() = default;

	virtual t_num eval(const ExprParser<t_num>& context) const override
	{
		const t_num val = m_child->eval(context);

		switch(m_op)
		{
			case '+': return val;
			case '-': return -val;
		}

		throw std::runtime_error{"Invalid unary operator."};
	}

	virtual void print(std::ostream& ostr=std::cout, std::size_t indent=0) const override
	{
		for(std::size_t i=0; i<indent; ++i) ostr << " | ";
		ostr << "unary operator " << m_op << "\n";
		m_child->print(ostr, indent+1);
	}

private:
	char m_op{'+'};
	std::shared_ptr<ExprAST<t_num>> m_child{};
};


template<typename t_num=double>
class ExprASTVar : public ExprAST<t_num>
{
public:
	ExprASTVar(const std::string& name) : m_name{name}
	{}

	ExprASTVar(const ExprASTVar<t_num>&) = delete;
	const ExprASTVar<t_num>& operator=(const ExprASTVar<t_num>&) = delete;

	virtual ~ExprASTVar() = default;

	virtual t_num eval(const ExprParser<t_num>& context) const override
	{
		return context.get_var(m_name);
	}

	virtual void print(std::ostream& ostr=std::cout, std::size_t indent=0) const override
	{
		for(std::size_t i=0; i<indent; ++i) ostr << " | ";
		ostr << "variable \"" << m_name << "\"\n";
	}

private:
	std::string m_name{};
};


template<typename t_num=double>
class ExprASTValue : public ExprAST<t_num>
{
public:
	ExprASTValue(t_num val) : m_val{val}
	{}

	ExprASTValue(const ExprASTValue<t_num>&) = delete;
	const ExprASTValue<t_num>& operator=(const ExprASTValue<t_num>&) = delete;

	virtual ~ExprASTValue() = default;

	virtual t_num eval(const ExprParser<t_num>&) const override
	{
		return m_val;
	}

	virtual void print(std::ostream& ostr=std::cout, std::size_t indent=0) const override
	{
		for(std::size_t i=0; i<indent; ++i) ostr << " | ";
		ostr << "value " << m_val << "\n";
	}

private:
	t_num m_val{};
};


template<typename t_num=double>
class ExprASTCall : public ExprAST<t_num>
{
public:
	ExprASTCall(const std::string& name)
		: m_name{name}
	{}

	ExprASTCall(const std::string& name, const std::shared_ptr<ExprAST<t_num>>& arg1)
		: m_name{name}
	{
		m_args.push_back(arg1);
	}

	ExprASTCall(const std::string& name,
		const std::shared_ptr<ExprAST<t_num>>& arg1, const std::shared_ptr<ExprAST<t_num>>& arg2)
		: m_name{name}
	{
		m_args.push_back(arg1);
		m_args.push_back(arg2);
	}

	ExprASTCall(const ExprASTCall<t_num>&) = delete;
	const ExprASTCall<t_num>& operator=(const ExprASTCall<t_num>&) = delete;

	virtual ~ExprASTCall() = default;

	virtual t_num eval(const ExprParser<t_num>& context) const override
	{
		if(m_args.size() == 0)
			return context.call_func0(m_name);
		else if(m_args.size() == 1)
			return context.call_func1(m_name, m_args[0]->eval(context));
		else if(m_args.size() == 2)
			return context.call_func2(m_name, m_args[0]->eval(context), m_args[1]->eval(context));

		throw std::runtime_error("Invalid function call.");
	}

	virtual void print(std::ostream& ostr=std::cout, std::size_t indent=0) const override
	{
		for(std::size_t i=0; i<indent; ++i) ostr << " | ";
		ostr << "function call \"" << m_name << "\"\n";

		for(const auto& arg : m_args)
			arg->print(ostr, indent+1);
	}

private:
	std::string m_name{};
	std::vector<std::shared_ptr<ExprAST<t_num>>> m_args{};
};


// ------------------------------------------------------------------------



template<typename t_num>
class ExprParser
{
	friend class ExprASTCall<t_num>;
	friend class ExprASTVar<t_num>;


public:
	ExprParser() : m_ast{}, m_vars{}, m_funcs0{}, m_funcs1{}, m_funcs2{}, m_istr{}, m_lookahead_text{}
	{
		register_funcs();
		register_consts();
	}


	bool parse(const std::string& str)
	{
		m_istr = std::make_shared<std::istringstream>(str);
		next_lookahead();
		m_ast = plus_term();

		//m_ast->print(std::cout);
		//std::cout << std::endl;

		// check if there would be are more tokens available?
		next_lookahead();
		bool at_eof = (m_lookahead == (int)Token::TOK_INVALID || m_lookahead == (int)Token::TOK_END);
		if(!at_eof)
			throw std::underflow_error("Not all input tokens have been consumed.");

		return !!m_ast;
	}


	t_num eval() const
	{
		if(!m_ast)
			throw std::runtime_error("Invalid AST.");
		return m_ast->eval(*this);
	}


protected:
	// ------------------------------------------------------------------------
	// tables / functions
	// ------------------------------------------------------------------------

	void register_funcs()
	{
		// common functions
		register_func1("abs", std::abs);
		register_func2("mod", modfunc<t_num>);

		// real functions
		if constexpr(std::is_floating_point_v<t_num>)
		{
			register_func1("sin", std::sin);
			register_func1("cos", std::cos);
			register_func1("tan", std::tan);
			register_func1("asin", std::asin);
			register_func1("acos", std::acos);
			register_func1("atan", std::atan);
			register_func1("sinh", std::sinh);
			register_func1("cosh", std::cosh);
			register_func1("tanh", std::tanh);
			register_func1("asinh", std::asinh);
			register_func1("acosh", std::acosh);
			register_func1("atanh", std::atanh);
			register_func1("sqrt", std::sqrt);
			register_func1("cbrt", std::cbrt);
			register_func1("exp", std::exp);
			register_func1("log", std::log);
			register_func1("log2", std::log2);
			register_func1("log10", std::log10);
			register_func1("erf", std::erf);
			register_func1("erfc", std::erfc);
			register_func1("round", std::round);
			register_func1("ceil", std::ceil);
			register_func1("floor", std::floor);
#ifdef HAS_ERFINV
			register_func1("erf_inv", boost::math::erf_inv);
#endif

			register_func2("pow", std::pow);
			register_func2("atan2", std::atan2);
		}

		// integer functions
		else if constexpr(std::is_integral_v<t_num>)
		{
			register_func2("pow", [](t_num t1, t_num t2) -> t_num { return t_num(std::pow(t1, t2)); } );
		}
	}


	// call function with zero parameters
	t_num call_func0(const std::string& strName) const
	{
		return m_funcs0.at(strName)();
	}


	// call function with one parameter
	t_num call_func1(const std::string& strName, t_num t) const
	{
		return m_funcs1.at(strName)(t);
	}


	// call function with two parameters
	t_num call_func2(const std::string& strName, t_num t1, t_num t2) const
	{
		return m_funcs2.at(strName)(t1, t2);
	}


	// register constants
	void register_consts()
	{
		// real constants
		if constexpr(std::is_floating_point_v<t_num>)
		{
			register_var("pi", __pi<t_num>);
			register_var("hbar",  t_num(hbar<t_num>/meV<t_num>/sec<t_num>));	// hbar in [meV s]
			register_var("kB",  t_num(kB<t_num>/meV<t_num>*kelvin<t_num>));	// kB in [meV / K]
		}

		// integer constants
		else if constexpr(std::is_integral_v<t_num>)
		{
		}
	}


	// get constant
	t_num get_var(const std::string& strName) const
	{
		return m_vars.at(strName);
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
	std::shared_ptr<ExprAST<t_num>> plus_term()
	{
		// plus_term -> mul_term plus_term_rest
		if(m_lookahead == '(' || m_lookahead == (int)Token::TOK_NUM || m_lookahead == (int)Token::TOK_IDENT)
		{
			auto term_val = mul_term();
			auto expr_rest_val = plus_term_rest(term_val);

			return expr_rest_val;
		}
		else if(m_lookahead == '+')	// unary +
		{
			next_lookahead();
			auto term_val = mul_term();
			auto expr_rest_val = plus_term_rest(term_val);

			return expr_rest_val;
		}
		else if(m_lookahead == '-')	// unary -
		{
			next_lookahead();
			auto term_val = std::make_shared<ExprASTUnOp<t_num>>('-', mul_term());
			auto expr_rest_val = plus_term_rest(term_val);

			return expr_rest_val;
		}

		std::ostringstream ostr;
		if(m_lookahead == 0 || m_lookahead == EOF)
			ostr << "EOF in " << __func__ << ".";
		else
			ostr << "Invalid lookahead in " << __func__ << ": " << m_lookahead << ".";
		throw std::runtime_error(ostr.str());
	}


	std::shared_ptr<ExprAST<t_num>> plus_term_rest(const std::shared_ptr<ExprAST<t_num>>& arg)
	{
		// plus_term_rest -> '+' mul_term plus_term_rest
		if(m_lookahead == '+')
		{
			next_lookahead();
			auto term_val = std::make_shared<ExprASTBinOp<t_num>>('+', arg, mul_term());
			auto expr_rest_val = plus_term_rest(term_val);

			return expr_rest_val;
		}

		// plus_term_rest -> '-' mul_term plus_term_rest
		else if(m_lookahead == '-')
		{
			next_lookahead();
			auto term_val = std::make_shared<ExprASTBinOp<t_num>>('-', arg, mul_term());
			auto expr_rest_val = plus_term_rest(term_val);

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
	}


	/**
	 * *,/,% terms
	 * (precedence 2)
	 */
	std::shared_ptr<ExprAST<t_num>> mul_term()
	{
		// mul_term -> pow_term mul_term_rest
		if(m_lookahead == '(' || m_lookahead == (int)Token::TOK_NUM || m_lookahead == (int)Token::TOK_IDENT)
		{
			auto factor_val = pow_term();
			auto term_rest_val = mul_term_rest(factor_val);

			return term_rest_val;
		}

		std::ostringstream ostr;
		ostr << "Invalid lookahead in " << __func__ << ": " << m_lookahead << ".";
		throw std::runtime_error(ostr.str());
	}


	std::shared_ptr<ExprAST<t_num>> mul_term_rest(const std::shared_ptr<ExprAST<t_num>>& arg)
	{
		// mul_term_rest -> '*' pow_term mul_term_rest
		if(m_lookahead == '*')
		{
			next_lookahead();
			auto factor_val = std::make_shared<ExprASTBinOp<t_num>>('*', arg, pow_term());
			auto term_rest_val = mul_term_rest(factor_val);

			return term_rest_val;
		}

		// mul_term_rest -> '/' pow_term mul_term_rest
		else if(m_lookahead == '/')
		{
			next_lookahead();
			auto factor_val = std::make_shared<ExprASTBinOp<t_num>>('/', arg, pow_term());
			auto term_rest_val = mul_term_rest(factor_val);

			return term_rest_val;
		}

		// mul_term_rest -> '%' pow_term mul_term_rest
		else if(m_lookahead == '%')
		{
			next_lookahead();
			auto factor_val = std::make_shared<ExprASTBinOp<t_num>>('%', arg, pow_term());
			auto term_rest_val = mul_term_rest(factor_val);

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
	}


	/**
	 * ^ terms
	 * (precedence 3)
	 */
	std::shared_ptr<ExprAST<t_num>> pow_term()
	{
		// pow_term -> factor pow_term_rest
		if(m_lookahead == '(' || m_lookahead == (int)Token::TOK_NUM || m_lookahead == (int)Token::TOK_IDENT)
		{
			auto factor_val = factor();
			auto term_rest_val = pow_term_rest(factor_val);

			return term_rest_val;
		}

		std::ostringstream ostr;
		ostr << "Invalid lookahead in " << __func__ << ": " << m_lookahead << ".";
		throw std::runtime_error(ostr.str());
	}


	std::shared_ptr<ExprAST<t_num>> pow_term_rest(const std::shared_ptr<ExprAST<t_num>>& arg)
	{
		// pow_term_rest -> '^' factor pow_term_rest
		if(m_lookahead == '^')
		{
			next_lookahead();
			auto factor_val = std::make_shared<ExprASTBinOp<t_num>>('^', arg, factor());
			auto term_rest_val = pow_term_rest(factor_val);

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
	}


	/**
	 * () terms, real factor or identifier
	 * (highest precedence, 4)
	 */
	std::shared_ptr<ExprAST<t_num>> factor()
	{
		// factor -> '(' plus_term ')'
		if(m_lookahead == '(')
		{
			next_lookahead();
			auto expr_val = plus_term();
			match(')');
			next_lookahead();

			return expr_val;
		}

		// factor -> TOK_NUM
		else if(m_lookahead == (int)Token::TOK_NUM)
		{
			t_num val = m_lookahead_val;
			next_lookahead();

			return std::make_shared<ExprASTValue<t_num>>(val);
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

					return std::make_shared<ExprASTCall<t_num>>(ident);
				}

				// function with arguments
				else
				{
					// first argument
					auto expr_val1 = plus_term();

					// one-argument-function
					// factor -> TOK_IDENT '(' plus_term ')'
					if(m_lookahead == ')')
					{
						next_lookahead();

						return std::make_shared<ExprASTCall<t_num>>(ident, expr_val1);
					}

					// two-argument-function
					// factor -> TOK_IDENT '(' plus_term ',' plus_term ')'
					else if(m_lookahead == ',')
					{
						next_lookahead();
						auto expr_val2 = plus_term();
						match(')');
						next_lookahead();

						return std::make_shared<ExprASTCall<t_num>>(ident, expr_val1, expr_val2);
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
				// register the variable if it doesn't yet exist
				if(m_vars.find(ident) == m_vars.end())
					register_var(ident, t_num{});
				return std::make_shared<ExprASTVar<t_num>>(ident);
			}
		}

		std::ostringstream ostr;
		ostr << "Invalid lookahead in " << __func__ << ": " << m_lookahead << ".";
		throw std::runtime_error(ostr.str());
	}
	// ----------------------------------------------------------------------------


public:
	// register a function with no parameters
	void register_func0(const std::string& name, t_num(*fkt)())
	{
		m_funcs0.emplace(std::make_pair(name, static_cast<t_num(*)()>(fkt)));
	}


	// register a function with one parameter
	void register_func1(const std::string& name, t_num(*fkt)(t_num))
	{
		m_funcs1.emplace(std::make_pair(name, static_cast<t_num(*)(t_num)>(fkt)));
	}


	// register a function with two parameters
	void register_func2(const std::string& name, t_num(*fkt)(t_num, t_num))
	{
		m_funcs2.emplace(std::make_pair(name, static_cast<t_num(*)(t_num, t_num)>(fkt)));
	}


	// register a variable or constant
	void register_var(const std::string& name, t_num val)
	{
		// overwrite value if key already exists
		if(auto [iter, ok] = m_vars.emplace(std::make_pair(name, val)); !ok)
			iter->second = val;
	}


	const std::unordered_map<std::string, t_num>& get_vars() const
	{
		return m_vars;
	}


private:
	// ast root
	std::shared_ptr<ExprAST<t_num>> m_ast{};

	// variables and constants
	std::unordered_map<std::string, t_num> m_vars{};

	// functions
	std::unordered_map<std::string, t_num(*)()> m_funcs0{};
	std::unordered_map<std::string, t_num(*)(t_num)> m_funcs1{};
	std::unordered_map<std::string, t_num(*)(t_num, t_num)> m_funcs2{};


private:
	std::shared_ptr<std::istream> m_istr{};

	int m_lookahead = (int)Token::TOK_INVALID;
	t_num m_lookahead_val = 0;
	std::string m_lookahead_text = "";
};


}
#endif
