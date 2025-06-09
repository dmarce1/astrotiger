#pragma once

#include <cassert>
#include <functional>
#include <memory>
#include <vector>

template<typename, typename >
struct Expression;

template<typename >
struct Valarray;

template<typename >
struct ValarrayExpression;

template<typename E>
concept Expr = requires(E e, int i) {
	typename E::value_type;      // must have a nested value_type
	{	e[i]}-> std::convertible_to<typename E::value_type>;
	{	e.size()}-> std::convertible_to<int>;
};

template<Expr Expression1, Expr Expression2, typename Operator, typename T>
struct BinaryOperator: public Expression<BinaryOperator<Expression1, Expression2, Operator, T>, T> {
	using value_type = T;
	BinaryOperator(Expression1 const&, Expression2 const&, Operator const&);
	T operator[](int) const;
	int size() const;
private:
	Expression1 expression1_{ };
	Expression2 expression2_{ };
	Operator operator_ { };
};

template<typename Derived, typename T>
struct Expression {
	using value_type = T;
	T operator[](int) const;
	int size() const;
};

template<typename T>
struct Valarray {
	using value_type = T;
	Valarray() = default;
	Valarray(Valarray<T> const &array);
	T operator[](int i) const;
	int size() const;
	operator ValarrayExpression<T>();
	friend class ValarrayExpression<T> ;
private:
	std::shared_ptr<std::vector<T>> dataPointer_;
};

template<typename T>
struct ValarrayExpression {
	using value_type = T;
	ValarrayExpression(Valarray<T> const &array);
	ValarrayExpression(Valarray<T> &array);
	T operator[](int i) const;
	int size() const;
private:
	std::shared_ptr<std::vector<T>> dataPointer_;
};

template<Expr Expression1, Expr Expression2, typename Operator, typename T>
BinaryOperator<Expression1, Expression2, Operator, T>::BinaryOperator(Expression1 const &expression1, Expression2 const &expression2, Operator const &_operator) :
		expression1_(expression1), expression2_(expression2), operator_(_operator) {
}

template<Expr Expression1, Expr Expression2, typename Operator, typename T>
BinaryOperator<Expression1, Expression2, Operator, T>::value_type BinaryOperator<Expression1, Expression2, Operator, T>::operator[](int i) const {
	return operator_(expression1_[i], expression2_[i]);
}

template<Expr Expression1, Expr Expression2, typename Operator, typename T>
int BinaryOperator<Expression1, Expression2, Operator, T>::size() const {
	return std::max(expression1_.size(), expression2_.size());
}

template<typename Derived, typename T>
T Expression<Derived, T>::operator[](int i) const {
	return static_cast<Derived const&>(*this)[i];
}

template<typename Derived, typename T>
int Expression<Derived, T>::size() const {
	return static_cast<Derived const&>(*this).size();
}

template<typename T>
Valarray<T>::Valarray(Valarray<T> const &array) :
		dataPointer_(array.dataPointer) {
}

template<typename T>
T Valarray<T>::operator[](int i) const {
	return (*dataPointer_)[i];
}

template<typename T>
int Valarray<T>::size() const {
	return dataPointer_->size();
}

template<typename T>
Valarray<T>::operator ValarrayExpression<T>() {
	return ValarrayExpression<T>(*this);
}

template<typename T>
ValarrayExpression<T>::ValarrayExpression(Valarray<T> const &array) :
		dataPointer_(array.dataPointer_) {
}

template<typename T>
T ValarrayExpression<T>::operator[](int i) const {
	return (*dataPointer_)[i];
}

template<typename T>
int ValarrayExpression<T>::size() const {
	return dataPointer_.size();
}

template<typename ...Args>
struct Valarray2Expression {
	template<typename T>
	T operator()(T const &expression) const {
		return expression;
	}
};

template<typename T>
struct Valarray2Expression<Valarray<T>> {
	auto operator()(Valarray<T> &expression) const {
		return ValarrayExpression<T>(expression);
	}
	auto operator()(Valarray<T> const &expression) const {
		return ValarrayExpression<T>(expression);
	}
};

template<Expr E1, Expr E2, typename Op>
BinaryOperator(E1 const&, E2 const&, Op const&) -> BinaryOperator<E1, E2, Op, typename E1::value_type>;

template<Expr E1, Expr E2>
auto operator+(E1 const &a, E2 const &b) {
	auto lhs = Valarray2Expression<E1>()(a);
	auto rhs = Valarray2Expression<E2>()(b);
	std::plus<typename E1::value_type> op;
	return BinaryOperator(lhs, rhs, op);  // CTAD deduces <decltype(lhs),decltype(rhs),…,value_type>
}

template<Expr Expression1, Expr Expression2>
auto operator-(Expression1 const &expression1, Expression2 const &expression2) {
	using T = typename Expression1::value_type;
	std::minus<T> _operator;
	return BinaryOperator<Expression1, Expression2, std::minus<T>, T>(expression1, expression2, _operator);
}

template<Expr Expression1, Expr Expression2>
auto operator*(Expression1 const &expression1, Expression2 const &expression2) {
	using T = typename Expression1::value_type;
	std::multiplies<T> _operator;
	return BinaryOperator<Expression1, Expression2, std::multiplies<T>, T>(expression1, expression2, _operator);
}

template<Expr Expression1, Expr Expression2>
auto operator/(Expression1 const &expression1, Expression2 const &expression2) {
	using T = typename Expression1::value_type;
	std::divides<T> _operator;
	return BinaryOperator<Expression1, Expression2, std::divides<T>, T>(expression1, expression2, _operator);
}

