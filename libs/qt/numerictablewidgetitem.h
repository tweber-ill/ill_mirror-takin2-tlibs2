/**
 * numeric table widget item
 * @author Tobias Weber <tweber@ill.fr>
 * @date Dec-2018
 * @license GPLv3, see 'LICENSE' file
 * @desc The present version was forked on 28-Dec-2018 and 21-Apr-2021 from my privately developed "misc" project (https://github.com/t-weber/misc).
 */

#ifndef __NUM_TABWIDGETITEM_H__
#define __NUM_TABWIDGETITEM_H__


#include <QtWidgets/QTableWidget>
#include "../str.h"


template<class T = double>
class NumericTableWidgetItem : public QTableWidgetItem
{
public:
	NumericTableWidgetItem(const T& val, std::streamsize prec = 6)
		: QTableWidgetItem(tl2::var_to_str(val, prec).c_str()),
		  m_val{val}, m_prec{prec}
	{}

	NumericTableWidgetItem(const QString& val, std::streamsize prec = 6)
		: QTableWidgetItem(val),
		  m_val{tl2::str_to_var<T>(val.toStdString())},
		  m_prec{prec}
	{}

	virtual bool operator<(const QTableWidgetItem& item) const override
	{
		T val1 = tl2::str_to_var<T>(text().toStdString());
		T val2 = tl2::str_to_var<T>(item.text().toStdString());

		return val1 < val2;
	}

	virtual QTableWidgetItem* clone() const override
	{
		auto item = new NumericTableWidgetItem<T>(m_val);
		item->setData(Qt::UserRole, this->data(Qt::UserRole));
		return item;
	};

	virtual void setData(int itemdatarole, const QVariant& var) override
	{
		m_val = var.value<T>();
		std::string str = tl2::var_to_str(m_val, m_prec);
		QTableWidgetItem::setData(itemdatarole, str.c_str());
	}

	T GetValue() const { return m_val; }

private:
	T m_val{};
	std::streamsize m_prec{6};
};


#endif
