/**
 * parser entry point
 * @author Tobias Weber <tweber@ill.fr>
 * @date 20-dec-19
 * @license see 'LICENSE' file
 * @desc Forked on 18/July/2020 from my privatly developed "matrix_calc" project (https://github.com/t-weber/matrix_calc).
 */

#include "ast.h"
#include "parser.h"
#include "llasm.h"
#include "printast.h"

#include <fstream>
#include <locale>
#include <boost/program_options.hpp>
namespace args = boost::program_options;


/**
 * Lexer message
 */
void yy::Lexer::LexerOutput(const char* str, int /*len*/)
{
	tl2::log_err("Lexer output (line ", GetCurLine(), "): ", str, ".");
}


/**
 * Lexer error output
 */
void yy::Lexer::LexerError(const char* err)
{
	tl2::log_err("Lexer error in line ", GetCurLine(), ": ", err, ".");
}


/**
 * Parser error output
 */
void yy::Parser::error(const std::string& err)
{
	tl2::log_err("Parser error in line ", context.GetCurLine(), ": ", err, ".");
}


/**
 * call lexer from parser
 */
extern yy::Parser::symbol_type yylex(yy::ParserContext &context)
{
	return context.GetLexer().lex();
}


int main(int argc, char** argv)
{
	try
	{
		std::ios_base::sync_with_stdio(0);
		std::locale loc{};
		std::locale::global(loc);

		tl2::log_info("--------------------------------------------------------------------------------");
		tl2::log_info("This is the tlibs2 scripting tool.");
		tl2::log_info("Author: Tobias Weber <tweber@ill.fr>, 2020.");
		tl2::log_info("Licensed under GPLv3.");
		tl2::log_info("--------------------------------------------------------------------------------");


		// llvm toolchain
		std::string tool_opt = "opt";
		std::string tool_bc = "llvm-as";
		std::string tool_bclink = "llvm-link";
		std::string tool_interp = "lli";
		std::string tool_s = "llc";
		std::string tool_o = "clang++";
		std::string tool_exec = "clang++";
		std::string tool_strip = "llvm-strip";


		// --------------------------------------------------------------------
		// get program arguments
		// --------------------------------------------------------------------
		std::vector<std::string> vecProgs;
		bool interpret = false;
		bool optimise = false;
		bool show_symbols = false;
		bool show_ast = false;
		std::string outprog;

		args::options_description arg_descr("Compiler arguments");
		arg_descr.add_options()
			("out,o", args::value(&outprog), "compiled program output")
			("optimise,O", args::bool_switch(&optimise), "optimise program")
			("interpret,i", args::bool_switch(&interpret), "directly run program in interpreter")
			("symbols,s", args::bool_switch(&show_symbols), "output symbol table")
			("ast,a", args::bool_switch(&show_ast), "output syntax tree")
			("program", args::value<decltype(vecProgs)>(&vecProgs), "input program to compile");

		args::positional_options_description posarg_descr;
		posarg_descr.add("program", -1);

		args::options_description arg_descr_toolchain("Toolchain programs");
		arg_descr_toolchain.add_options()
			("tool_opt", args::value(&tool_opt), "llvm optimiser")
			("tool_bc", args::value(&tool_bc), "llvm bitcode assembler")
			("tool_bclink", args::value(&tool_bclink), "llvm bitcode linker")
			("tool_interp", args::value(&tool_interp), "llvm bitcode interpreter")
			("tool_bccomp", args::value(&tool_s), "llvm bitcode compiler")
			("tool_asm", args::value(&tool_o), "native assembler")
			("tool_link", args::value(&tool_exec), "native linker")
			("tool_strip", args::value(&tool_strip), "strip tool");
		arg_descr.add(arg_descr_toolchain);

		auto argparser = args::command_line_parser{argc, argv};
		argparser.style(args::command_line_style::default_style);
		argparser.options(arg_descr);
		argparser.positional(posarg_descr);

		args::variables_map mapArgs;
		auto parsedArgs = argparser.run();
		args::store(parsedArgs, mapArgs);
		args::notify(mapArgs);

		if(vecProgs.size() == 0)
		{
			tl2::log_info("Please specify an input program.");
			tl2::log_info(arg_descr);
			return 0;
		}

		if(outprog == "")
		{
			outprog = "out";
			tl2::log_warn("No program output specified, using \"", outprog, "\".");
		}

		std::string outprog_ast = outprog + "_ast.xml";
		std::string outprog_syms = outprog + "_syms.txt";

		std::string outprog_3ac = outprog + ".asm";
		std::string outprog_3ac_opt = outprog + "_opt.asm";
		std::string outprog_bc = outprog + ".bc";
		std::string outprog_linkedbc = outprog + "_linked.bc";
		std::string outprog_s = outprog + ".s";
		std::string outprog_o = outprog + ".o";

		std::string runtime_3ac = optimise ? "runtime_opt.asm" : "runtime.asm";
		std::string runtime_bc = "runtime.bc";
		// --------------------------------------------------------------------



		// --------------------------------------------------------------------
		// parse input
		// --------------------------------------------------------------------
		const std::string& inprog = vecProgs[0];
		tl2::log_info("Parsing \"", inprog, "\"...");

		std::ifstream ifstr{inprog};
		if(!ifstr)
		{
			tl2::log_err("Cannot open \"", inprog, "\".");
			return -1;
		}
		yy::ParserContext ctx{ifstr};

		// register runtime functions
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "pow", SymbolType::SCALAR, {SymbolType::SCALAR, SymbolType::SCALAR});
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "sin", SymbolType::SCALAR, {SymbolType::SCALAR});
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "cos", SymbolType::SCALAR, {SymbolType::SCALAR});
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "sqrt", SymbolType::SCALAR, {SymbolType::SCALAR});
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "exp", SymbolType::SCALAR, {SymbolType::SCALAR});
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "fabs", SymbolType::SCALAR, {SymbolType::SCALAR});
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "labs", SymbolType::INT, {SymbolType::INT});

		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "set_eps", SymbolType::VOID, {SymbolType::SCALAR});
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "get_eps", SymbolType::SCALAR, {});

		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "strlen", SymbolType::INT, {SymbolType::STRING});

		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "putstr", SymbolType::VOID, {SymbolType::STRING});
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "putflt", SymbolType::VOID, {SymbolType::SCALAR});
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "putint", SymbolType::VOID, {SymbolType::INT});
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "getflt", SymbolType::SCALAR, {SymbolType::STRING});
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "getint", SymbolType::INT, {SymbolType::STRING});

		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "flt_to_str", SymbolType::VOID, {SymbolType::SCALAR, SymbolType::STRING, SymbolType::INT});
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "int_to_str", SymbolType::VOID, {SymbolType::INT, SymbolType::STRING, SymbolType::INT});


		yy::Parser parser(ctx);
		int res = parser.parse();
		if(res != 0)
		{
			tl2::log_err("Parser reports failure.");
			return res;
		}

		if(show_symbols)
		{
			tl2::log_info("Writing symbol table to \"", outprog_syms, "\"...");

			std::ofstream ostrSyms{outprog_syms};
			ostrSyms << ctx.GetSymbols() << std::endl;
		}

		if(show_ast)
		{
			tl2::log_info("Writing AST to \"", outprog_ast, "\"...");

			std::ofstream ostrAST{outprog_ast};
			ASTPrinter printer{&ostrAST};

			ostrAST << "<ast>\n";
			auto stmts = ctx.GetStatements()->GetStatementList();
			for(auto iter=stmts.rbegin(); iter!=stmts.rend(); ++iter)
			{
				(*iter)->accept(&printer);
				ostrAST << "\n";
			}
			ostrAST << "</ast>" << std::endl;
		}
		// --------------------------------------------------------------------



		// --------------------------------------------------------------------
		// 3AC generation
		// --------------------------------------------------------------------
		tl2::log_info("Generating intermediate code: \"",
			inprog, "\" -> \"", outprog_3ac, "\"...");

		std::ofstream ofstr{outprog_3ac};
		std::ostream* ostr = &ofstr;
		LLAsm llasm{&ctx.GetSymbols(), ostr};
		auto stmts = ctx.GetStatements()->GetStatementList();
		for(auto iter=stmts.rbegin(); iter!=stmts.rend(); ++iter)
		{
			(*iter)->accept(&llasm);
			(*ostr) << std::endl;
		}


		// additional runtime/startup code
		(*ostr) << "\n" << R"START(
; -----------------------------------------------------------------------------
; imported libc functions
declare double @pow(double, double)
declare double @sin(double)
declare double @cos(double)
declare double @sqrt(double)
declare double @exp(double)
declare double @fabs(double)
declare i64 @labs(i64)

declare i64 @strlen(i8*)
declare i8* @strncpy(i8*, i8*, i64)
declare i8* @strncat(i8*, i8*, i64)
declare i32 @strncmp(i8*, i8*, i64)
declare i32 @puts(i8*)
declare i32 @snprintf(i8*, i64, i8*, ...)
declare i32 @printf(i8*, ...)
declare i32 @scanf(i8*, ...)
declare i8* @memcpy(i8*, i8*, i64)
declare i8* @malloc(i64)
declare i8* @calloc(i64, i64)
declare void @free(i8*)
; -----------------------------------------------------------------------------


; -----------------------------------------------------------------------------
; external runtime functions from runtime.c
declare void @ext_set_eps(double)
declare double @ext_get_eps()

declare double @ext_determinant(double*, i64)
declare i64 @ext_power(double*, double*, i64, i64)
declare i64 @ext_transpose(double*, double*, i64, i64)
; -----------------------------------------------------------------------------


; -----------------------------------------------------------------------------
; constants
@__strfmt_s = constant [3 x i8] c"%s\00"
@__strfmt_lg = constant [4 x i8] c"%lg\00"
@__strfmt_ld = constant [4 x i8] c"%ld\00"
@__str_vecbegin = constant [3 x i8] c"[ \00"
@__str_vecend = constant [3 x i8] c" ]\00"
@__str_vecsep = constant [3 x i8] c", \00"
@__str_matsep = constant [3 x i8] c"; \00"
; -----------------------------------------------------------------------------


; -----------------------------------------------------------------------------
; runtime functions

; get the user epsilon
define double @get_eps()
{
	%eps = call double @ext_get_eps()
	ret double %eps
}

; set the user epsilon
define void @set_eps(double %eps)
{
	call void (double) @ext_set_eps(double %eps)
	ret void
}

; returns 0 if flt <= eps
define double @zero_eps(double %flt)
{
	%eps = call double @get_eps()
	%fltabs = call double (double) @fabs(double %flt)

	%cond = fcmp ole double %fltabs, %eps
	br i1 %cond, label %labelIf, label %labelEnd
labelIf:
	ret double 0.
labelEnd:
	ret double %flt
}

; double -> string
define void @flt_to_str(double %flt, i8* %strptr, i64 %len)
{
	%fmtptr = bitcast [4 x i8]* @__strfmt_lg to i8*
	%theflt = call double (double) @zero_eps(double %flt)
	call i32 (i8*, i64, i8*, ...) @snprintf(i8* %strptr, i64 %len, i8* %fmtptr, double %theflt)
	ret void
}

; int -> string
define void @int_to_str(i64 %i, i8* %strptr, i64 %len)
{
	%fmtptr = bitcast [4 x i8]* @__strfmt_ld to i8*
	call i32 (i8*, i64, i8*, ...) @snprintf(i8* %strptr, i64 %len, i8* %fmtptr, i64 %i)
	ret void
}

; output a string
define void @putstr(i8* %val)
{
	call i32 (i8*) @puts(i8* %val)
	ret void
}

; output a float
define void @putflt(double %val)
{
	; convert to string
	%strval = alloca [64 x i8]
	%strvalptr = bitcast [64 x i8]* %strval to i8*
	call void @flt_to_str(double %val, i8* %strvalptr, i64 64)

	; output string
	call void (i8*) @putstr(i8* %strvalptr)
	ret void
}

; output an int
define void @putint(i64 %val)
{
	; convert to string
	%strval = alloca [64 x i8]
	%strvalptr = bitcast [64 x i8]* %strval to i8*
	call void @int_to_str(i64 %val, i8* %strvalptr, i64 64)

	; output string
	call void (i8*) @putstr(i8* %strvalptr)
	ret void
}

; input a float
define double @getflt(i8* %str)
{
	; output given string
	%fmtptr_s = bitcast [3 x i8]* @__strfmt_s to i8*
	call i32 (i8*, ...) @printf(i8* %fmtptr_s, i8* %str)

	; alloc double
	%d_ptr = alloca double

	; read double from stdin
	%fmtptr_g = bitcast [4 x i8]* @__strfmt_lg to i8*
	call i32 (i8*, ...) @scanf(i8* %fmtptr_g, double* %d_ptr)

	%d = load double, double* %d_ptr
	ret double %d
}

; input an int
define i64 @getint(i8* %str)
{
	; output given string
	%fmtptr_s = bitcast [3 x i8]* @__strfmt_s to i8*
	call i32 (i8*, ...) @printf(i8* %fmtptr_s, i8* %str)

	; alloc int
	%i_ptr = alloca i64

	; read int from stdin
	%fmtptr_ld = bitcast [4 x i8]* @__strfmt_ld to i8*
	call i32 (i8*, ...) @scanf(i8* %fmtptr_ld, i64* %i_ptr)

	%i = load i64, i64* %i_ptr
	ret i64 %i
}

; -----------------------------------------------------------------------------


; -----------------------------------------------------------------------------
; main entry point for llvm
define i32 @main()
{
	; call entry function
	call void @start()

	ret i32 0
}
; -----------------------------------------------------------------------------
)START";

		(*ostr) << std::endl;
		// --------------------------------------------------------------------



		// --------------------------------------------------------------------
		// 3AC optimisation
		// --------------------------------------------------------------------
		if(optimise)
		{
			tl2::log_info("Optimising intermediate code: \"",
				outprog_3ac, "\" -> \"", outprog_3ac_opt, "\"...");

			std::string cmd_opt = tool_opt + " -stats -S --strip-debug -o "
				+ outprog_3ac_opt + " " + outprog_3ac;
			if(std::system(cmd_opt.c_str()) != 0)
			{
				tl2::log_err("Failed.");
				return -1;
			}

			outprog_3ac = outprog_3ac_opt;
		}
		// --------------------------------------------------------------------



		// --------------------------------------------------------------------
		// Bitcode generation
		// --------------------------------------------------------------------
		tl2::log_info("Assembling bitcode: \"",
			outprog_3ac, "\" -> \"", outprog_bc, "\"...");

		std::string cmd_bc = tool_bc + " -o " + outprog_bc + " " + outprog_3ac;
		if(std::system(cmd_bc.c_str()) != 0)
		{
			tl2::log_err("Failed.");
			return -1;
		}


		tl2::log_info("Assembling runtime bitcode: \"",
			runtime_3ac, "\" -> \"", runtime_bc, "\"...");

		cmd_bc = tool_bc + " -o " + runtime_bc + " " + runtime_3ac;
		if(std::system(cmd_bc.c_str()) != 0)
		{
			tl2::log_err("Failed.");
			return -1;
		}
		// --------------------------------------------------------------------



		// --------------------------------------------------------------------
		// Bitcode linking
		// --------------------------------------------------------------------
		tl2::log_info("Linking bitcode to runtime: \"",
			outprog_bc, "\" + \"", runtime_bc, "\" -> \"",
			outprog_linkedbc, "\"...");

		std::string cmd_bclink = tool_bclink + " -o " + outprog_linkedbc + " " + outprog_bc + " " + runtime_bc;
		if(std::system(cmd_bclink.c_str()) != 0)
		{
			tl2::log_err("Failed.");
			return -1;
		}
		// --------------------------------------------------------------------


		// interpret bitcode
		if(interpret)
		{
			tl2::log_info("Interpreting bitcode \"", outprog_linkedbc, "\"...");

			std::string cmd_interp = tool_interp + " " + outprog_linkedbc;
			if(std::system(cmd_interp.c_str()) != 0)
			{
				tl2::log_err("Failed.");
				return -1;
			}
		}

		// compile bitcode
		else
		{
			tl2::log_info("Generating native assembly \"",
				outprog_linkedbc, "\" -> \"", outprog_s, "\"...");

			std::string opt_flag_s = optimise ? "-O2" : "";
			std::string cmd_s = tool_s + " " + opt_flag_s + " -o " + outprog_s + " " + outprog_linkedbc;
			if(std::system(cmd_s.c_str()) != 0)
			{
				tl2::log_err("Failed.");
				return -1;
			}


			tl2::log_info("Assembling native code \"", outprog_s, "\" -> \"", outprog_o, "\"...");

			std::string opt_flag_o = optimise ? "-O2" : "";
			std::string cmd_o = tool_o + " " + opt_flag_o + " -c -o " + outprog_o + " " + outprog_s;
			if(std::system(cmd_o.c_str()) != 0)
			{
				tl2::log_err("Failed.");
				return -1;
			}


			tl2::log_info("Generating native executable \"",
				outprog_o, "\" -> \"", outprog, "\"...");

			std::string opt_flag_exec = optimise ? "-O2" : "";
			std::string cmd_exec = tool_exec + " " + opt_flag_exec + " -o " + outprog + " " + outprog_o;
			if(std::system(cmd_exec.c_str()) != 0)
			{
				tl2::log_err("Failed.");
				return -1;
			}


			if(optimise)
			{
				tl2::log_info("Stripping debug symbols from \"", outprog, "\"...");

				std::string cmd_strip = tool_strip + " " + outprog;
				if(std::system(cmd_strip.c_str()) != 0)
				{
					tl2::log_err("Failed.");
					return -1;
				}
			}
		}
	}
	catch(const std::exception& ex)
	{
		tl2::log_err("Error: ", ex.what());
		return -1;
	}

	return 0;
}
