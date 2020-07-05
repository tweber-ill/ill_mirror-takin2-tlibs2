/**
 * llvm three-address code generator
 * @author Tobias Weber <tweber@ill.fr>
 * @date 11-apr-20
 * @license see 'LICENSE' file
 * @desc Forked on 5-July-2020 from the privately developed "matrix_calc" project (https://github.com/t-weber/matrix_calc).
 *
 * References:
 *	* https://llvm.org/docs/tutorial/MyFirstLanguageFrontend/LangImpl03.html
 *	* https://llvm.org/docs/GettingStarted.html
 *	* https://llvm.org/docs/LangRef.html
 */

#ifndef __LLASM_H__
#define __LLASM_H__

#include "ast.h"
#include "sym.h"


class LLAsm : public ASTVisitor
{
public:
	LLAsm(SymTab* syms, std::ostream* ostr=&std::cout);
	virtual ~LLAsm() = default;


	virtual t_astret visit(const ASTUMinus* ast) override;
	virtual t_astret visit(const ASTPlus* ast) override;
	virtual t_astret visit(const ASTMult* ast) override;
	virtual t_astret visit(const ASTMod* ast) override;
	virtual t_astret visit(const ASTPow* ast) override;
	virtual t_astret visit(const ASTTransp* ast) override;
	virtual t_astret visit(const ASTNorm* ast) override;
	virtual t_astret visit(const ASTVar* ast) override;
	virtual t_astret visit(const ASTCall* ast) override;
	virtual t_astret visit(const ASTStmts* ast) override;
	virtual t_astret visit(const ASTVarDecl* ast) override;
	virtual t_astret visit(const ASTFunc* ast) override;
	virtual t_astret visit(const ASTReturn* ast) override;
	virtual t_astret visit(const ASTAssign* ast) override;
	virtual t_astret visit(const ASTArrayAccess* ast) override;
	virtual t_astret visit(const ASTArrayAssign* ast) override;
	virtual t_astret visit(const ASTComp* ast) override;
	virtual t_astret visit(const ASTCond* ast) override;
	virtual t_astret visit(const ASTBool* ast) override;
	virtual t_astret visit(const ASTLoop* ast) override;
	virtual t_astret visit(const ASTStrConst* ast) override;
	virtual t_astret visit(const ASTExprList* ast) override;
	virtual t_astret visit(const ASTNumConst<double>* ast) override;
	virtual t_astret visit(const ASTNumConst<std::int64_t>* ast) override;


	// ------------------------------------------------------------------------
	// internally handled dummy nodes
	// ------------------------------------------------------------------------
	virtual t_astret visit(const ASTArgNames*) override { return nullptr; }
	virtual t_astret visit(const ASTTypeDecl*) override { return nullptr; }
	// ------------------------------------------------------------------------


protected:
	t_astret get_tmp_var(SymbolType ty = SymbolType::SCALAR,
		const std::array<std::size_t, 2>* dims = nullptr,
		const std::string* name = nullptr, bool on_heap = false);

	/**
	 * create a label for branch instructions
	 */
	std::string get_label();

	/**
	 * get the type name for a symbol
	 */
	static std::string get_type_name(SymbolType ty);

	/**
	 * find the symbol with a specific name in the symbol table
	 */
	t_astret get_sym(const std::string& name) const;

	/**
	 * convert symbol to another type
	 */
	t_astret convert_sym(t_astret sym, SymbolType ty_to);


	/**
	 * generates if-then-else code
	 */
	template<class t_funcCond, class t_funcBody, class t_funcElseBody>
	void generate_cond(t_funcCond funcCond, t_funcBody funcBody, t_funcElseBody funcElseBody, bool hasElse=0);


	/**
	 * generates loop code
	 */
	template<class t_funcCond, class t_funcBody>
	void generate_loop(t_funcCond funcCond, t_funcBody funcBody);


	/**
	 * generates loop code with managed counter
	 */
	template<class t_funcBody>
	void generate_loop(std::int64_t start, std::int64_t end, t_funcBody funcBody);


private:
	std::size_t m_varCount = 0;	// # of tmp vars
	std::size_t m_labelCount = 0;	// # of labels

	std::vector<std::string> m_curscope;

	SymTab* m_syms = nullptr;

	std::ostream* m_ostr = &std::cout;

	// helper functions to reduce code redundancy
	t_astret scalar_matrix_prod(t_astret scalar, t_astret matrix, bool mul_or_div=1);
};



/**
 * generates if-then-else code
 */
template<class t_funcCond, class t_funcBody, class t_funcElseBody>
void LLAsm::generate_cond(t_funcCond funcCond, t_funcBody funcBody, t_funcElseBody funcElseBody, bool hasElse)
{
	(*m_ostr) << "\n;-------------------------------------------------------------\n";
	(*m_ostr) << "; condition head\n";
	(*m_ostr) << ";-------------------------------------------------------------\n";
	t_astret cond = funcCond();
	(*m_ostr) << ";-------------------------------------------------------------\n";

	std::string labelIf = get_label();
	std::string labelElse = hasElse ? get_label() : "";
	std::string labelEnd = get_label();

	if(hasElse)
		(*m_ostr) << "br i1 %" << cond->name << ", label %" << labelIf << ", label %" << labelElse << "\n";
	else
		(*m_ostr) << "br i1 %" << cond->name << ", label %" << labelIf << ", label %" << labelEnd << "\n";

	(*m_ostr) << ";-------------------------------------------------------------\n";
	(*m_ostr) << "; condition body\n";
	(*m_ostr) << ";-------------------------------------------------------------\n";
	(*m_ostr) << labelIf << ":\n";
	funcBody();
	(*m_ostr) << ";-------------------------------------------------------------\n";

	(*m_ostr) << "br label %" << labelEnd << "\n";

	if(hasElse)
	{
		(*m_ostr) << ";-------------------------------------------------------------\n";
		(*m_ostr) << "; condition \"else\" body\n";
		(*m_ostr) << ";-------------------------------------------------------------\n";
		(*m_ostr) << labelElse << ":\n";
		funcElseBody();
		(*m_ostr) << ";-------------------------------------------------------------\n";
		(*m_ostr) << "br label %" << labelEnd << "\n";
	}

	(*m_ostr) << labelEnd << ":\n";
	(*m_ostr) << ";-------------------------------------------------------------\n\n";
}



/**
 * generates loop code
 */
template<class t_funcCond, class t_funcBody>
void LLAsm::generate_loop(t_funcCond funcCond, t_funcBody funcBody)
{
	std::string labelStart = get_label();
	std::string labelBegin = get_label();
	std::string labelEnd = get_label();

	(*m_ostr) << "\n;-------------------------------------------------------------\n";
	(*m_ostr) << "; loop head\n";
	(*m_ostr) << ";-------------------------------------------------------------\n";
	(*m_ostr) << "br label %" << labelStart << "\n";
	(*m_ostr) << labelStart << ":\n";
	t_astret cond = funcCond();
	(*m_ostr) << "br i1 %" << cond->name << ", label %" << labelBegin << ", label %" << labelEnd << "\n";

	(*m_ostr) << ";-------------------------------------------------------------\n";
	(*m_ostr) << "; loop body\n";
	(*m_ostr) << ";-------------------------------------------------------------\n";
	(*m_ostr) << labelBegin << ":\n";
	funcBody();
	(*m_ostr) << ";-------------------------------------------------------------\n";

	(*m_ostr) << "br label %" << labelStart << "\n";
	(*m_ostr) << labelEnd << ":\n";
	(*m_ostr) << ";-------------------------------------------------------------\n\n";
}



/**
 * generates loop code with managed counter
 */
template<class t_funcBody>
void LLAsm::generate_loop(std::int64_t start, std::int64_t end, t_funcBody funcBody)
{
	(*m_ostr) << "\n;-------------------------------------------------------------\n";
	(*m_ostr) << "; loop counter\n";
	(*m_ostr) << ";-------------------------------------------------------------\n";
	t_astret ctr = get_tmp_var(SymbolType::INT);
	t_astret ctrval = get_tmp_var(SymbolType::INT);
	(*m_ostr) << "%" << ctr->name << " = alloca i64\n";
	(*m_ostr) << "store i64 " << start << ", i64* %" << ctr->name << "\n";

	generate_loop([this, ctr, ctrval, end]() -> t_astret
	{
		// loop condition: ctr < end
		(*m_ostr) << "%" << ctrval->name << " = load i64, i64* %" << ctr->name << "\n";

		t_astret _cond = get_tmp_var();
		(*m_ostr) << "%" << _cond->name << " = icmp slt i64 %" << ctrval->name <<  ", " << end << "\n";

		return _cond;
	}, [this, &funcBody, ctr, ctrval]
	{
		funcBody(ctrval);

		(*m_ostr) << ";-------------------------------------------------------------\n";
		(*m_ostr) << "; increment loop counter\n";
		(*m_ostr) << ";-------------------------------------------------------------\n";
		t_astret newctrval = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << newctrval->name << " = add i64 %" << ctrval->name << ", 1\n";
		(*m_ostr) << "store i64 %" << newctrval->name << ", i64* %" << ctr->name << "\n";
	});
}

#endif
