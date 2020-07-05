/**
 * lexer
 * @author Tobias Weber <tweber@ill.fr>
 * @date 09-dec-19
 * @license see 'LICENSE' file
 * @desc Forked on 5-July-2020 from the privately developed "matrix_calc" project (https://github.com/t-weber/matrix_calc).
 *
 * References:
 *	https://github.com/westes/flex/tree/master/examples/manual
 *	http://www.gnu.org/software/bison/manual/html_node/index.html
 *	http://git.savannah.gnu.org/cgit/bison.git/tree/examples
 */

%option c++
%option noyywrap
%option yyclass = "yy::Lexer"


%{
	#include <string>
	#include "parser.h"
%}


white		[ \t\r]
sign		[+-]
num		[0-9]
int		{num}+
sint		{sign}?{int}
real		{int}\.{num}*|\.{int}
ident		[A-Za-z_]+[A-Za-z0-9_]*


%%

"\n"		{ IncCurLine(); }
{white}+	/* whitespace */
"#".*		/* comment */

";"			{ return yytext[0]; }
"+"|"-"		{ return yytext[0]; }
"*"|"/"|"%"	{ return yytext[0]; }
"^"|"'"			{ return yytext[0]; }
"("|")"		{ return yytext[0]; }
"{"|"}"		{ return yytext[0]; }
"["|"]"		{ return yytext[0]; }
"|"|"|"		{ return yytext[0]; }
","			{ return yytext[0]; }

"=="		{ return yy::Parser::make_EQU(); }
"<>"		{ return yy::Parser::make_NEQ(); }
">="		{ return yy::Parser::make_GEQ(); }
"<="		{ return yy::Parser::make_LEQ(); }
">"			{ return yy::Parser::make_GT(); }
"<"			{ return yy::Parser::make_LT(); }
"and"		{ return yy::Parser::make_AND(); }
"xor"		{ return yy::Parser::make_XOR(); }
"or"		{ return yy::Parser::make_OR(); }
"not"		{ return yy::Parser::make_NOT(); }

"="			{ return yytext[0]; }

"scalar"|"var"	{ return yy::Parser::make_SCALARDECL(); }
"vec"		{ return yy::Parser::make_VECTORDECL(); }
"mat"		{ return yy::Parser::make_MATRIXDECL(); }
"str"		{ return yy::Parser::make_STRINGDECL(); }
"int"		{ return yy::Parser::make_INTDECL(); }

"func"		{ return yy::Parser::make_FUNC(); }
"ret"		{ return yy::Parser::make_RET(); }

"if"		{ return yy::Parser::make_IF(); }
"then"		{ return yy::Parser::make_THEN(); }
"else"		{ return yy::Parser::make_ELSE(); }
"loop"		{ return yy::Parser::make_LOOP(); }
"do"		{ return yy::Parser::make_DO(); }

{int}		{ return yy::Parser::make_INT(std::stol(yytext)); }
{real}([eE]{sint})?	{ return yy::Parser::make_REAL(std::stod(yytext)); }
{ident}		{ return yy::Parser::make_IDENT(yytext); }
"\""[^\"]*"\""	{
	std::string str{yytext+1, yytext+yyleng-1};
	return yy::Parser::make_STRING(str); }

. {
	LexerError(
		(std::string{"Invalid token: \""} + 
		std::string{yytext} + 
		std::string{"\"."}).c_str());
}
%%
