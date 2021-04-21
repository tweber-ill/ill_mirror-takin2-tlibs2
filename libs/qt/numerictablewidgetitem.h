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
	NumericTableWidgetItem(T&& val)
		: QTableWidgetItem(std::to_string(std::forward<T>(val)).c_str()),
		  m_val{val}
	{}

	NumericTableWidgetItem(const T& val)
		: QTableWidgetItem(std::to_string(val).c_str()),
		  m_val{val}
	{}

	NumericTableWidgetItem(const QString& val)
		: QTableWidgetItem(val),
		  m_val{tl2::str_to_var<T>(val.toStdString())}
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
		QTableWidgetItem::setData(itemdatarole, std::to_string(m_val).c_str());
	}

	T GetValue() const { return m_val; }

private:
	T m_val{};
};


#endif
