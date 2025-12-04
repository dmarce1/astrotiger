#include <array>


template<char C>
struct FreeIndex {
   static constexpr char value = C;
};


template<typename, int>
struct TensorExpression0;

template<typename, int, char>
struct TensorExpression1;

template<typename, int, char, char>
struct TensorExpression2;

template<typename, int, char, char, char>
struct TensorExpression3;

template<typename, int, char, char, char, char>
struct TensorExpression4;


template<typename H, int D>
struct TensorExpression0 {
   TensorExpression0(H const& h) : handle_(h) {
   }
   template<typename F>
   auto& operator=(TensorExpression0<F, D> const& other) {
      (*this)() = other();
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression0<F, D> const& other) {
      (*this)() += other();
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression0<F, D> const& other) {
      (*this)() -= other();
      return *this;
   }
private:
   H handle_;
};


template<typename H, int D, char I>
struct TensorExpression1 {
   TensorExpression1(H const& h) : handle_(h) {
   }
   template<typename F>
   auto& operator=(TensorExpression1<F, D, I> const& other) {
      for(int i; i < D; i++) {
         (*this)(i) = other(i);
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression1<F, D, I> const& other) {
      for(int i; i < D; i++) {
         (*this)(i) += other(i);
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression1<F, D, I> const& other) {
      for(int i; i < D; i++) {
         (*this)(i) -= other(i);
      }
      return *this;
   }
   template<typename F, char J>
   auto operator*(TensorExpression1<F, D, J> const& other);
   template<typename F>
   auto operator*(TensorExpression1<F, D, I> const& other);
   template<typename F, char J, char K>
   auto operator*(TensorExpression2<F, D, J, K> const& other);
   template<typename F, char J>
   auto operator*(TensorExpression2<F, D, I, J> const& other);
   template<typename F, char J>
   auto operator*(TensorExpression2<F, D, J, I> const& other);
   template<typename F, char J, char K, char L>
   auto operator*(TensorExpression3<F, D, J, K, L> const& other);
   template<typename F, char J, char K>
   auto operator*(TensorExpression3<F, D, I, J, K> const& other);
   template<typename F, char J, char K>
   auto operator*(TensorExpression3<F, D, J, I, K> const& other);
   template<typename F, char J, char K>
   auto operator*(TensorExpression3<F, D, J, K, I> const& other);
   template<typename F, char J, char K, char L>
   auto operator*(TensorExpression4<F, D, I, J, K, L> const& other);
   template<typename F, char J, char K, char L>
   auto operator*(TensorExpression4<F, D, J, I, K, L> const& other);
   template<typename F, char J, char K, char L>
   auto operator*(TensorExpression4<F, D, J, K, I, L> const& other);
   template<typename F, char J, char K, char L>
   auto operator*(TensorExpression4<F, D, J, K, L, I> const& other);
private:
   H handle_;
};


template<typename H, int D, char I, char J>
struct TensorExpression2 {
   TensorExpression2(H const& h) : handle_(h) {
   }
   template<typename F>
   auto& operator=(TensorExpression2<F, D, I, J> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            (*this)(i, j) = other(i, j);
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator=(TensorExpression2<F, D, J, I> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            (*this)(i, j) = other(j, i);
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression2<F, D, I, J> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            (*this)(i, j) += other(i, j);
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression2<F, D, J, I> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            (*this)(i, j) += other(j, i);
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression2<F, D, I, J> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            (*this)(i, j) -= other(i, j);
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression2<F, D, J, I> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            (*this)(i, j) -= other(j, i);
         }
      }
      return *this;
   }
   template<typename F, char K>
   auto operator*(TensorExpression1<F, D, K> const& other);
   template<typename F>
   auto operator*(TensorExpression1<F, D, I> const& other);
   template<typename F>
   auto operator*(TensorExpression1<F, D, J> const& other);
   template<typename F, char K, char L>
   auto operator*(TensorExpression2<F, D, K, L> const& other);
   template<typename F, char K>
   auto operator*(TensorExpression2<F, D, I, K> const& other);
   template<typename F, char K>
   auto operator*(TensorExpression2<F, D, K, I> const& other);
   template<typename F, char K>
   auto operator*(TensorExpression2<F, D, J, K> const& other);
   template<typename F, char K>
   auto operator*(TensorExpression2<F, D, K, J> const& other);
   template<typename F>
   auto operator*(TensorExpression2<F, D, I, J> const& other);
   template<typename F>
   auto operator*(TensorExpression2<F, D, J, I> const& other);
   template<typename F, char K, char L>
   auto operator*(TensorExpression3<F, D, I, K, L> const& other);
   template<typename F, char K, char L>
   auto operator*(TensorExpression3<F, D, K, I, L> const& other);
   template<typename F, char K, char L>
   auto operator*(TensorExpression3<F, D, K, L, I> const& other);
   template<typename F, char K, char L>
   auto operator*(TensorExpression3<F, D, J, K, L> const& other);
   template<typename F, char K, char L>
   auto operator*(TensorExpression3<F, D, K, J, L> const& other);
   template<typename F, char K, char L>
   auto operator*(TensorExpression3<F, D, K, L, J> const& other);
   template<typename F, char K>
   auto operator*(TensorExpression3<F, D, I, J, K> const& other);
   template<typename F, char K>
   auto operator*(TensorExpression3<F, D, I, K, J> const& other);
   template<typename F, char K>
   auto operator*(TensorExpression3<F, D, J, I, K> const& other);
   template<typename F, char K>
   auto operator*(TensorExpression3<F, D, J, K, I> const& other);
   template<typename F, char K>
   auto operator*(TensorExpression3<F, D, K, I, J> const& other);
   template<typename F, char K>
   auto operator*(TensorExpression3<F, D, K, J, I> const& other);
   template<typename F, char K, char L, char M>
   auto operator*(TensorExpression4<F, D, I, K, L, M> const& other);
   template<typename F, char K, char L, char M>
   auto operator*(TensorExpression4<F, D, K, I, L, M> const& other);
   template<typename F, char K, char L, char M>
   auto operator*(TensorExpression4<F, D, K, L, I, M> const& other);
   template<typename F, char K, char L, char M>
   auto operator*(TensorExpression4<F, D, K, L, M, I> const& other);
   template<typename F, char K, char L, char M>
   auto operator*(TensorExpression4<F, D, J, K, L, M> const& other);
   template<typename F, char K, char L, char M>
   auto operator*(TensorExpression4<F, D, K, J, L, M> const& other);
   template<typename F, char K, char L, char M>
   auto operator*(TensorExpression4<F, D, K, L, J, M> const& other);
   template<typename F, char K, char L, char M>
   auto operator*(TensorExpression4<F, D, K, L, M, J> const& other);
   template<typename F, char K, char L>
   auto operator*(TensorExpression4<F, D, I, J, K, L> const& other);
   template<typename F, char K, char L>
   auto operator*(TensorExpression4<F, D, I, K, J, L> const& other);
   template<typename F, char K, char L>
   auto operator*(TensorExpression4<F, D, I, K, L, J> const& other);
   template<typename F, char K, char L>
   auto operator*(TensorExpression4<F, D, J, I, K, L> const& other);
   template<typename F, char K, char L>
   auto operator*(TensorExpression4<F, D, J, K, I, L> const& other);
   template<typename F, char K, char L>
   auto operator*(TensorExpression4<F, D, J, K, L, I> const& other);
   template<typename F, char K, char L>
   auto operator*(TensorExpression4<F, D, K, I, J, L> const& other);
   template<typename F, char K, char L>
   auto operator*(TensorExpression4<F, D, K, I, L, J> const& other);
   template<typename F, char K, char L>
   auto operator*(TensorExpression4<F, D, K, J, I, L> const& other);
   template<typename F, char K, char L>
   auto operator*(TensorExpression4<F, D, K, J, L, I> const& other);
   template<typename F, char K, char L>
   auto operator*(TensorExpression4<F, D, K, L, I, J> const& other);
   template<typename F, char K, char L>
   auto operator*(TensorExpression4<F, D, K, L, J, I> const& other);
private:
   H handle_;
};


template<typename H, int D, char I, char J, char K>
struct TensorExpression3 {
   TensorExpression3(H const& h) : handle_(h) {
   }
   template<typename F>
   auto& operator=(TensorExpression3<F, D, I, J, K> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               (*this)(i, j, k) = other(i, j, k);
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator=(TensorExpression3<F, D, I, K, J> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               (*this)(i, j, k) = other(i, k, j);
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator=(TensorExpression3<F, D, J, I, K> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               (*this)(i, j, k) = other(j, i, k);
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator=(TensorExpression3<F, D, J, K, I> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               (*this)(i, j, k) = other(j, k, i);
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator=(TensorExpression3<F, D, K, I, J> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               (*this)(i, j, k) = other(k, i, j);
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator=(TensorExpression3<F, D, K, J, I> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               (*this)(i, j, k) = other(k, j, i);
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression3<F, D, I, J, K> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               (*this)(i, j, k) += other(i, j, k);
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression3<F, D, I, K, J> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               (*this)(i, j, k) += other(i, k, j);
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression3<F, D, J, I, K> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               (*this)(i, j, k) += other(j, i, k);
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression3<F, D, J, K, I> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               (*this)(i, j, k) += other(j, k, i);
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression3<F, D, K, I, J> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               (*this)(i, j, k) += other(k, i, j);
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression3<F, D, K, J, I> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               (*this)(i, j, k) += other(k, j, i);
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression3<F, D, I, J, K> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               (*this)(i, j, k) -= other(i, j, k);
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression3<F, D, I, K, J> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               (*this)(i, j, k) -= other(i, k, j);
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression3<F, D, J, I, K> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               (*this)(i, j, k) -= other(j, i, k);
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression3<F, D, J, K, I> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               (*this)(i, j, k) -= other(j, k, i);
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression3<F, D, K, I, J> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               (*this)(i, j, k) -= other(k, i, j);
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression3<F, D, K, J, I> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               (*this)(i, j, k) -= other(k, j, i);
            }
         }
      }
      return *this;
   }
   template<typename F, char L>
   auto operator*(TensorExpression1<F, D, L> const& other);
   template<typename F>
   auto operator*(TensorExpression1<F, D, I> const& other);
   template<typename F>
   auto operator*(TensorExpression1<F, D, J> const& other);
   template<typename F>
   auto operator*(TensorExpression1<F, D, K> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression2<F, D, I, L> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression2<F, D, L, I> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression2<F, D, J, L> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression2<F, D, L, J> const& other);
   template<typename F>
   auto operator*(TensorExpression2<F, D, I, J> const& other);
   template<typename F>
   auto operator*(TensorExpression2<F, D, J, I> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression2<F, D, K, L> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression2<F, D, L, K> const& other);
   template<typename F>
   auto operator*(TensorExpression2<F, D, I, K> const& other);
   template<typename F>
   auto operator*(TensorExpression2<F, D, K, I> const& other);
   template<typename F>
   auto operator*(TensorExpression2<F, D, J, K> const& other);
   template<typename F>
   auto operator*(TensorExpression2<F, D, K, J> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression3<F, D, I, L, M> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression3<F, D, L, I, M> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression3<F, D, L, M, I> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression3<F, D, J, L, M> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression3<F, D, L, J, M> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression3<F, D, L, M, J> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression3<F, D, I, J, L> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression3<F, D, I, L, J> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression3<F, D, J, I, L> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression3<F, D, J, L, I> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression3<F, D, L, I, J> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression3<F, D, L, J, I> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression3<F, D, K, L, M> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression3<F, D, L, K, M> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression3<F, D, L, M, K> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression3<F, D, I, K, L> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression3<F, D, I, L, K> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression3<F, D, K, I, L> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression3<F, D, K, L, I> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression3<F, D, L, I, K> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression3<F, D, L, K, I> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression3<F, D, J, K, L> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression3<F, D, J, L, K> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression3<F, D, K, J, L> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression3<F, D, K, L, J> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression3<F, D, L, J, K> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression3<F, D, L, K, J> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, I, J, K> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, I, K, J> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, J, I, K> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, J, K, I> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, K, I, J> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, K, J, I> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, I, J, L, M> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, I, L, J, M> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, I, L, M, J> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, J, I, L, M> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, J, L, I, M> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, J, L, M, I> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, L, I, J, M> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, L, I, M, J> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, L, J, I, M> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, L, J, M, I> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, L, M, I, J> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, L, M, J, I> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, I, K, L, M> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, I, L, K, M> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, I, L, M, K> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, K, I, L, M> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, K, L, I, M> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, K, L, M, I> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, L, I, K, M> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, L, I, M, K> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, L, K, I, M> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, L, K, M, I> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, L, M, I, K> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, L, M, K, I> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, J, K, L, M> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, J, L, K, M> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, J, L, M, K> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, K, J, L, M> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, K, L, J, M> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, K, L, M, J> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, L, J, K, M> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, L, J, M, K> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, L, K, J, M> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, L, K, M, J> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, L, M, J, K> const& other);
   template<typename F, char L, char M>
   auto operator*(TensorExpression4<F, D, L, M, K, J> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression4<F, D, I, J, K, L> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression4<F, D, I, J, L, K> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression4<F, D, I, K, J, L> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression4<F, D, I, K, L, J> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression4<F, D, I, L, J, K> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression4<F, D, I, L, K, J> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression4<F, D, J, I, K, L> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression4<F, D, J, I, L, K> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression4<F, D, J, K, I, L> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression4<F, D, J, K, L, I> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression4<F, D, J, L, I, K> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression4<F, D, J, L, K, I> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression4<F, D, K, I, J, L> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression4<F, D, K, I, L, J> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression4<F, D, K, J, I, L> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression4<F, D, K, J, L, I> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression4<F, D, K, L, I, J> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression4<F, D, K, L, J, I> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression4<F, D, L, I, J, K> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression4<F, D, L, I, K, J> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression4<F, D, L, J, I, K> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression4<F, D, L, J, K, I> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression4<F, D, L, K, I, J> const& other);
   template<typename F, char L>
   auto operator*(TensorExpression4<F, D, L, K, J, I> const& other);
private:
   H handle_;
};


template<typename H, int D, char I, char J, char K, char L>
struct TensorExpression4 {
   TensorExpression4(H const& h) : handle_(h) {
   }
   template<typename F>
   auto& operator=(TensorExpression4<F, D, I, J, K, L> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) = other(i, j, k, l);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator=(TensorExpression4<F, D, I, J, L, K> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) = other(i, j, l, k);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator=(TensorExpression4<F, D, I, K, J, L> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) = other(i, k, j, l);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator=(TensorExpression4<F, D, I, K, L, J> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) = other(i, k, l, j);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator=(TensorExpression4<F, D, I, L, J, K> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) = other(i, l, j, k);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator=(TensorExpression4<F, D, I, L, K, J> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) = other(i, l, k, j);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator=(TensorExpression4<F, D, J, I, K, L> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) = other(j, i, k, l);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator=(TensorExpression4<F, D, J, I, L, K> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) = other(j, i, l, k);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator=(TensorExpression4<F, D, J, K, I, L> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) = other(j, k, i, l);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator=(TensorExpression4<F, D, J, K, L, I> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) = other(j, k, l, i);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator=(TensorExpression4<F, D, J, L, I, K> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) = other(j, l, i, k);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator=(TensorExpression4<F, D, J, L, K, I> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) = other(j, l, k, i);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator=(TensorExpression4<F, D, K, I, J, L> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) = other(k, i, j, l);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator=(TensorExpression4<F, D, K, I, L, J> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) = other(k, i, l, j);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator=(TensorExpression4<F, D, K, J, I, L> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) = other(k, j, i, l);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator=(TensorExpression4<F, D, K, J, L, I> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) = other(k, j, l, i);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator=(TensorExpression4<F, D, K, L, I, J> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) = other(k, l, i, j);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator=(TensorExpression4<F, D, K, L, J, I> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) = other(k, l, j, i);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator=(TensorExpression4<F, D, L, I, J, K> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) = other(l, i, j, k);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator=(TensorExpression4<F, D, L, I, K, J> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) = other(l, i, k, j);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator=(TensorExpression4<F, D, L, J, I, K> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) = other(l, j, i, k);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator=(TensorExpression4<F, D, L, J, K, I> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) = other(l, j, k, i);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator=(TensorExpression4<F, D, L, K, I, J> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) = other(l, k, i, j);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator=(TensorExpression4<F, D, L, K, J, I> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) = other(l, k, j, i);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression4<F, D, I, J, K, L> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) += other(i, j, k, l);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression4<F, D, I, J, L, K> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) += other(i, j, l, k);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression4<F, D, I, K, J, L> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) += other(i, k, j, l);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression4<F, D, I, K, L, J> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) += other(i, k, l, j);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression4<F, D, I, L, J, K> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) += other(i, l, j, k);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression4<F, D, I, L, K, J> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) += other(i, l, k, j);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression4<F, D, J, I, K, L> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) += other(j, i, k, l);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression4<F, D, J, I, L, K> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) += other(j, i, l, k);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression4<F, D, J, K, I, L> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) += other(j, k, i, l);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression4<F, D, J, K, L, I> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) += other(j, k, l, i);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression4<F, D, J, L, I, K> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) += other(j, l, i, k);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression4<F, D, J, L, K, I> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) += other(j, l, k, i);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression4<F, D, K, I, J, L> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) += other(k, i, j, l);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression4<F, D, K, I, L, J> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) += other(k, i, l, j);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression4<F, D, K, J, I, L> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) += other(k, j, i, l);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression4<F, D, K, J, L, I> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) += other(k, j, l, i);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression4<F, D, K, L, I, J> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) += other(k, l, i, j);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression4<F, D, K, L, J, I> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) += other(k, l, j, i);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression4<F, D, L, I, J, K> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) += other(l, i, j, k);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression4<F, D, L, I, K, J> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) += other(l, i, k, j);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression4<F, D, L, J, I, K> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) += other(l, j, i, k);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression4<F, D, L, J, K, I> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) += other(l, j, k, i);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression4<F, D, L, K, I, J> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) += other(l, k, i, j);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator+=(TensorExpression4<F, D, L, K, J, I> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) += other(l, k, j, i);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression4<F, D, I, J, K, L> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) -= other(i, j, k, l);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression4<F, D, I, J, L, K> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) -= other(i, j, l, k);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression4<F, D, I, K, J, L> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) -= other(i, k, j, l);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression4<F, D, I, K, L, J> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) -= other(i, k, l, j);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression4<F, D, I, L, J, K> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) -= other(i, l, j, k);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression4<F, D, I, L, K, J> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) -= other(i, l, k, j);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression4<F, D, J, I, K, L> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) -= other(j, i, k, l);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression4<F, D, J, I, L, K> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) -= other(j, i, l, k);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression4<F, D, J, K, I, L> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) -= other(j, k, i, l);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression4<F, D, J, K, L, I> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) -= other(j, k, l, i);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression4<F, D, J, L, I, K> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) -= other(j, l, i, k);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression4<F, D, J, L, K, I> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) -= other(j, l, k, i);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression4<F, D, K, I, J, L> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) -= other(k, i, j, l);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression4<F, D, K, I, L, J> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) -= other(k, i, l, j);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression4<F, D, K, J, I, L> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) -= other(k, j, i, l);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression4<F, D, K, J, L, I> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) -= other(k, j, l, i);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression4<F, D, K, L, I, J> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) -= other(k, l, i, j);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression4<F, D, K, L, J, I> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) -= other(k, l, j, i);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression4<F, D, L, I, J, K> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) -= other(l, i, j, k);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression4<F, D, L, I, K, J> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) -= other(l, i, k, j);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression4<F, D, L, J, I, K> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) -= other(l, j, i, k);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression4<F, D, L, J, K, I> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) -= other(l, j, k, i);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression4<F, D, L, K, I, J> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) -= other(l, k, i, j);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto& operator-=(TensorExpression4<F, D, L, K, J, I> const& other) {
      for(int i; i < D; i++) {
         for(int j; j < D; j++) {
            for(int k; k < D; k++) {
               for(int l; l < D; l++) {
                  (*this)(i, j, k, l) -= other(l, k, j, i);
               }
            }
         }
      }
      return *this;
   }
   template<typename F>
   auto operator*(TensorExpression1<F, D, I> const& other);
   template<typename F>
   auto operator*(TensorExpression1<F, D, J> const& other);
   template<typename F>
   auto operator*(TensorExpression1<F, D, K> const& other);
   template<typename F>
   auto operator*(TensorExpression1<F, D, L> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression2<F, D, I, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression2<F, D, M, I> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression2<F, D, J, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression2<F, D, M, J> const& other);
   template<typename F>
   auto operator*(TensorExpression2<F, D, I, J> const& other);
   template<typename F>
   auto operator*(TensorExpression2<F, D, J, I> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression2<F, D, K, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression2<F, D, M, K> const& other);
   template<typename F>
   auto operator*(TensorExpression2<F, D, I, K> const& other);
   template<typename F>
   auto operator*(TensorExpression2<F, D, K, I> const& other);
   template<typename F>
   auto operator*(TensorExpression2<F, D, J, K> const& other);
   template<typename F>
   auto operator*(TensorExpression2<F, D, K, J> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression2<F, D, L, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression2<F, D, M, L> const& other);
   template<typename F>
   auto operator*(TensorExpression2<F, D, I, L> const& other);
   template<typename F>
   auto operator*(TensorExpression2<F, D, L, I> const& other);
   template<typename F>
   auto operator*(TensorExpression2<F, D, J, L> const& other);
   template<typename F>
   auto operator*(TensorExpression2<F, D, L, J> const& other);
   template<typename F>
   auto operator*(TensorExpression2<F, D, K, L> const& other);
   template<typename F>
   auto operator*(TensorExpression2<F, D, L, K> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, I, J, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, I, M, J> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, J, I, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, J, M, I> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, M, I, J> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, M, J, I> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, I, K, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, I, M, K> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, K, I, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, K, M, I> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, M, I, K> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, M, K, I> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, J, K, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, J, M, K> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, K, J, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, K, M, J> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, M, J, K> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, M, K, J> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, I, J, K> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, I, K, J> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, J, I, K> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, J, K, I> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, K, I, J> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, K, J, I> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, I, L, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, I, M, L> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, L, I, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, L, M, I> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, M, I, L> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, M, L, I> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, J, L, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, J, M, L> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, L, J, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, L, M, J> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, M, J, L> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, M, L, J> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, I, J, L> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, I, L, J> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, J, I, L> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, J, L, I> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, L, I, J> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, L, J, I> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, K, L, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, K, M, L> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, L, K, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, L, M, K> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, M, K, L> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression3<F, D, M, L, K> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, I, K, L> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, I, L, K> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, K, I, L> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, K, L, I> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, L, I, K> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, L, K, I> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, J, K, L> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, J, L, K> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, K, J, L> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, K, L, J> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, L, J, K> const& other);
   template<typename F>
   auto operator*(TensorExpression3<F, D, L, K, J> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, I, J, M, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, I, M, J, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, I, M, N, J> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, J, I, M, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, J, M, I, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, J, M, N, I> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, I, J, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, I, N, J> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, J, I, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, J, N, I> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, N, I, J> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, N, J, I> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, I, K, M, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, I, M, K, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, I, M, N, K> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, K, I, M, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, K, M, I, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, K, M, N, I> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, I, K, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, I, N, K> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, K, I, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, K, N, I> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, N, I, K> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, N, K, I> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, J, K, M, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, J, M, K, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, J, M, N, K> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, K, J, M, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, K, M, J, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, K, M, N, J> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, J, K, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, J, N, K> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, K, J, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, K, N, J> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, N, J, K> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, N, K, J> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, I, J, K, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, I, J, M, K> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, I, K, J, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, I, K, M, J> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, I, M, J, K> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, I, M, K, J> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, J, I, K, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, J, I, M, K> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, J, K, I, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, J, K, M, I> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, J, M, I, K> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, J, M, K, I> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, K, I, J, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, K, I, M, J> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, K, J, I, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, K, J, M, I> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, K, M, I, J> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, K, M, J, I> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, M, I, J, K> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, M, I, K, J> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, M, J, I, K> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, M, J, K, I> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, M, K, I, J> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, M, K, J, I> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, I, L, M, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, I, M, L, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, I, M, N, L> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, L, I, M, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, L, M, I, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, L, M, N, I> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, I, L, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, I, N, L> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, L, I, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, L, N, I> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, N, I, L> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, N, L, I> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, J, L, M, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, J, M, L, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, J, M, N, L> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, L, J, M, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, L, M, J, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, L, M, N, J> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, J, L, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, J, N, L> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, L, J, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, L, N, J> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, N, J, L> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, N, L, J> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, I, J, L, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, I, J, M, L> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, I, L, J, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, I, L, M, J> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, I, M, J, L> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, I, M, L, J> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, J, I, L, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, J, I, M, L> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, J, L, I, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, J, L, M, I> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, J, M, I, L> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, J, M, L, I> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, L, I, J, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, L, I, M, J> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, L, J, I, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, L, J, M, I> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, L, M, I, J> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, L, M, J, I> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, M, I, J, L> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, M, I, L, J> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, M, J, I, L> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, M, J, L, I> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, M, L, I, J> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, M, L, J, I> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, K, L, M, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, K, M, L, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, K, M, N, L> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, L, K, M, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, L, M, K, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, L, M, N, K> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, K, L, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, K, N, L> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, L, K, N> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, L, N, K> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, N, K, L> const& other);
   template<typename F, char M, char N>
   auto operator*(TensorExpression4<F, D, M, N, L, K> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, I, K, L, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, I, K, M, L> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, I, L, K, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, I, L, M, K> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, I, M, K, L> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, I, M, L, K> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, K, I, L, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, K, I, M, L> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, K, L, I, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, K, L, M, I> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, K, M, I, L> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, K, M, L, I> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, L, I, K, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, L, I, M, K> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, L, K, I, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, L, K, M, I> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, L, M, I, K> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, L, M, K, I> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, M, I, K, L> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, M, I, L, K> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, M, K, I, L> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, M, K, L, I> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, M, L, I, K> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, M, L, K, I> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, J, K, L, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, J, K, M, L> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, J, L, K, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, J, L, M, K> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, J, M, K, L> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, J, M, L, K> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, K, J, L, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, K, J, M, L> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, K, L, J, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, K, L, M, J> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, K, M, J, L> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, K, M, L, J> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, L, J, K, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, L, J, M, K> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, L, K, J, M> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, L, K, M, J> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, L, M, J, K> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, L, M, K, J> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, M, J, K, L> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, M, J, L, K> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, M, K, J, L> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, M, K, L, J> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, M, L, J, K> const& other);
   template<typename F, char M>
   auto operator*(TensorExpression4<F, D, M, L, K, J> const& other);
   template<typename F>
   auto operator*(TensorExpression4<F, D, I, J, K, L> const& other);
   template<typename F>
   auto operator*(TensorExpression4<F, D, I, J, L, K> const& other);
   template<typename F>
   auto operator*(TensorExpression4<F, D, I, K, J, L> const& other);
   template<typename F>
   auto operator*(TensorExpression4<F, D, I, K, L, J> const& other);
   template<typename F>
   auto operator*(TensorExpression4<F, D, I, L, J, K> const& other);
   template<typename F>
   auto operator*(TensorExpression4<F, D, I, L, K, J> const& other);
   template<typename F>
   auto operator*(TensorExpression4<F, D, J, I, K, L> const& other);
   template<typename F>
   auto operator*(TensorExpression4<F, D, J, I, L, K> const& other);
   template<typename F>
   auto operator*(TensorExpression4<F, D, J, K, I, L> const& other);
   template<typename F>
   auto operator*(TensorExpression4<F, D, J, K, L, I> const& other);
   template<typename F>
   auto operator*(TensorExpression4<F, D, J, L, I, K> const& other);
   template<typename F>
   auto operator*(TensorExpression4<F, D, J, L, K, I> const& other);
   template<typename F>
   auto operator*(TensorExpression4<F, D, K, I, J, L> const& other);
   template<typename F>
   auto operator*(TensorExpression4<F, D, K, I, L, J> const& other);
   template<typename F>
   auto operator*(TensorExpression4<F, D, K, J, I, L> const& other);
   template<typename F>
   auto operator*(TensorExpression4<F, D, K, J, L, I> const& other);
   template<typename F>
   auto operator*(TensorExpression4<F, D, K, L, I, J> const& other);
   template<typename F>
   auto operator*(TensorExpression4<F, D, K, L, J, I> const& other);
   template<typename F>
   auto operator*(TensorExpression4<F, D, L, I, J, K> const& other);
   template<typename F>
   auto operator*(TensorExpression4<F, D, L, I, K, J> const& other);
   template<typename F>
   auto operator*(TensorExpression4<F, D, L, J, I, K> const& other);
   template<typename F>
   auto operator*(TensorExpression4<F, D, L, J, K, I> const& other);
   template<typename F>
   auto operator*(TensorExpression4<F, D, L, K, I, J> const& other);
   template<typename F>
   auto operator*(TensorExpression4<F, D, L, K, J, I> const& other);
private:
   H handle_;
};


template<typename T, int D>
struct Tensor1 {
   T operator()(int i) const {
      return data_[flatten(i)];
   }
   T& operator()(int i) {
      return data_[flatten(i)];
   }
   template<char I>
   auto operator()(FreeIndex<I>) {
      return Tensor1Expression<Tensor1&, D, I>(*this);
   }
private:
   static constexpr int size = D;
   int flatten(int i) const {
      return i;
   }
   std::array<T, size> data_;
};


template<typename T, int D>
struct Tensor2 {
   T operator()(int i, int j) const {
      return data_[flatten(i, j)];
   }
   T& operator()(int i, int j) {
      return data_[flatten(i, j)];
   }
   template<char I, char J>
   auto operator()(FreeIndex<I>, FreeIndex<J>) {
      return Tensor2Expression<Tensor2&, D, I, J>(*this);
   }
   template<char I>
   auto operator()(FreeIndex<I>, FreeIndex<I>) {
      auto const lambda = [this]() {
         T sum = T(0);
         for( int i; i < D; i++ ) {
            sum += (*this)(i, i);
         }
         return sum;
      };
      return Tensor0Expression<decltype(lambda), D>(*this);
   }
private:
   static constexpr int size = D * D;
   int flatten(int i, int j) const {
      return j + D * i;
   }
   std::array<T, size> data_;
};


template<typename T, int D>
struct Tensor3 {
   T operator()(int i, int j, int k) const {
      return data_[flatten(i, j, k)];
   }
   T& operator()(int i, int j, int k) {
      return data_[flatten(i, j, k)];
   }
   template<char I, char J, char K>
   auto operator()(FreeIndex<I>, FreeIndex<J>, FreeIndex<K>) {
      return Tensor3Expression<Tensor3&, D, I, J, K>(*this);
   }
   template<char I, char J>
   auto operator()(FreeIndex<I>, FreeIndex<I>, FreeIndex<J>) {
      auto const lambda = [this](int j) {
         T sum = T(0);
         for( int i; i < D; i++ ) {
            sum += (*this)(i, i, j);
         }
         return sum;
      };
      return Tensor1Expression<decltype(lambda), D, J>(*this);
   }
   template<char I, char J>
   auto operator()(FreeIndex<I>, FreeIndex<J>, FreeIndex<I>) {
      auto const lambda = [this](int j) {
         T sum = T(0);
         for( int i; i < D; i++ ) {
            sum += (*this)(i, j, i);
         }
         return sum;
      };
      return Tensor1Expression<decltype(lambda), D, J>(*this);
   }
   template<char I, char J>
   auto operator()(FreeIndex<J>, FreeIndex<I>, FreeIndex<I>) {
      auto const lambda = [this](int j) {
         T sum = T(0);
         for( int i; i < D; i++ ) {
            sum += (*this)(j, i, i);
         }
         return sum;
      };
      return Tensor1Expression<decltype(lambda), D, J>(*this);
   }
private:
   static constexpr int size = D * D * D;
   int flatten(int i, int j, int k) const {
      return k + D * (j + D * i);
   }
   std::array<T, size> data_;
};


template<typename T, int D>
struct Tensor4 {
   T operator()(int i, int j, int k, int l) const {
      return data_[flatten(i, j, k, l)];
   }
   T& operator()(int i, int j, int k, int l) {
      return data_[flatten(i, j, k, l)];
   }
   template<char I, char J, char K, char L>
   auto operator()(FreeIndex<I>, FreeIndex<J>, FreeIndex<K>, FreeIndex<L>) {
      return Tensor4Expression<Tensor4&, D, I, J, K, L>(*this);
   }
   template<char I, char J, char K>
   auto operator()(FreeIndex<I>, FreeIndex<I>, FreeIndex<J>, FreeIndex<K>) {
      auto const lambda = [this](int j, int k) {
         T sum = T(0);
         for( int i; i < D; i++ ) {
            sum += (*this)(i, i, j, k);
         }
         return sum;
      };
      return Tensor2Expression<decltype(lambda), D, J, K>(*this);
   }
   template<char I, char J, char K>
   auto operator()(FreeIndex<I>, FreeIndex<J>, FreeIndex<I>, FreeIndex<K>) {
      auto const lambda = [this](int j, int k) {
         T sum = T(0);
         for( int i; i < D; i++ ) {
            sum += (*this)(i, j, i, k);
         }
         return sum;
      };
      return Tensor2Expression<decltype(lambda), D, J, K>(*this);
   }
   template<char I, char J, char K>
   auto operator()(FreeIndex<I>, FreeIndex<J>, FreeIndex<K>, FreeIndex<I>) {
      auto const lambda = [this](int j, int k) {
         T sum = T(0);
         for( int i; i < D; i++ ) {
            sum += (*this)(i, j, k, i);
         }
         return sum;
      };
      return Tensor2Expression<decltype(lambda), D, J, K>(*this);
   }
   template<char I, char J, char K>
   auto operator()(FreeIndex<J>, FreeIndex<I>, FreeIndex<I>, FreeIndex<K>) {
      auto const lambda = [this](int j, int k) {
         T sum = T(0);
         for( int i; i < D; i++ ) {
            sum += (*this)(j, i, i, k);
         }
         return sum;
      };
      return Tensor2Expression<decltype(lambda), D, J, K>(*this);
   }
   template<char I, char J, char K>
   auto operator()(FreeIndex<J>, FreeIndex<I>, FreeIndex<K>, FreeIndex<I>) {
      auto const lambda = [this](int j, int k) {
         T sum = T(0);
         for( int i; i < D; i++ ) {
            sum += (*this)(j, i, k, i);
         }
         return sum;
      };
      return Tensor2Expression<decltype(lambda), D, J, K>(*this);
   }
   template<char I, char J, char K>
   auto operator()(FreeIndex<J>, FreeIndex<K>, FreeIndex<I>, FreeIndex<I>) {
      auto const lambda = [this](int j, int k) {
         T sum = T(0);
         for( int i; i < D; i++ ) {
            sum += (*this)(j, k, i, i);
         }
         return sum;
      };
      return Tensor2Expression<decltype(lambda), D, J, K>(*this);
   }
   template<char I, char J>
   auto operator()(FreeIndex<I>, FreeIndex<I>, FreeIndex<J>, FreeIndex<J>) {
      auto const lambda = [this]() {
         T sum = T(0);
         for( int i; i < D; i++ ) {
            for( int j; j < D; j++ ) {
               sum += (*this)(i, i, j, j);
            }
         }
         return sum;
      };
      return Tensor0Expression<decltype(lambda), D>(*this);
   }
   template<char I, char J>
   auto operator()(FreeIndex<I>, FreeIndex<J>, FreeIndex<I>, FreeIndex<J>) {
      auto const lambda = [this]() {
         T sum = T(0);
         for( int i; i < D; i++ ) {
            for( int j; j < D; j++ ) {
               sum += (*this)(i, j, i, j);
            }
         }
         return sum;
      };
      return Tensor0Expression<decltype(lambda), D>(*this);
   }
   template<char I, char J>
   auto operator()(FreeIndex<I>, FreeIndex<J>, FreeIndex<J>, FreeIndex<I>) {
      auto const lambda = [this]() {
         T sum = T(0);
         for( int i; i < D; i++ ) {
            for( int j; j < D; j++ ) {
               sum += (*this)(i, j, j, i);
            }
         }
         return sum;
      };
      return Tensor0Expression<decltype(lambda), D>(*this);
   }
private:
   static constexpr int size = D * D * D * D;
   int flatten(int i, int j, int k, int l) const {
      return l + D * (k + D * (j + D * i));
   }
   std::array<T, size> data_;
};


template<typename H, int D, char I>
template<typename F, char J>
auto TensorExpression1<H, D, I>::operator*(TensorExpression1<F, D, J> const& other) {
   auto const lambda = [this, other](int i, int j) {
      using type = std::remove_cvref_t<decltype((*this)(0) * other(0))>;
      auto sum = type(0);
      sum += (*this)(i) * other(j);
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, J>(lambda);
}

template<typename H, int D, char I>
template<typename F>
auto TensorExpression1<H, D, I>::operator*(TensorExpression1<F, D, I> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0) * other(0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         sum += (*this)(i) * other(i);
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I>
template<typename F, char J, char K>
auto TensorExpression1<H, D, I>::operator*(TensorExpression2<F, D, J, K> const& other) {
   auto const lambda = [this, other](int i, int j, int k) {
      using type = std::remove_cvref_t<decltype((*this)(0) * other(0, 0))>;
      auto sum = type(0);
      sum += (*this)(i) * other(j, k);
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, J, K>(lambda);
}

template<typename H, int D, char I>
template<typename F, char J>
auto TensorExpression1<H, D, I>::operator*(TensorExpression2<F, D, I, J> const& other) {
   auto const lambda = [this, other](int j) {
      using type = std::remove_cvref_t<decltype((*this)(0) * other(0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         sum += (*this)(i) * other(i, j);
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, J>(lambda);
}

template<typename H, int D, char I>
template<typename F, char J>
auto TensorExpression1<H, D, I>::operator*(TensorExpression2<F, D, J, I> const& other) {
   auto const lambda = [this, other](int j) {
      using type = std::remove_cvref_t<decltype((*this)(0) * other(0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         sum += (*this)(i) * other(j, i);
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, J>(lambda);
}

template<typename H, int D, char I>
template<typename F, char J, char K, char L>
auto TensorExpression1<H, D, I>::operator*(TensorExpression3<F, D, J, K, L> const& other) {
   auto const lambda = [this, other](int i, int j, int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0) * other(0, 0, 0))>;
      auto sum = type(0);
      sum += (*this)(i) * other(j, k, l);
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, J, K, L>(lambda);
}

template<typename H, int D, char I>
template<typename F, char J, char K>
auto TensorExpression1<H, D, I>::operator*(TensorExpression3<F, D, I, J, K> const& other) {
   auto const lambda = [this, other](int j, int k) {
      using type = std::remove_cvref_t<decltype((*this)(0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         sum += (*this)(i) * other(i, j, k);
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, K>(lambda);
}

template<typename H, int D, char I>
template<typename F, char J, char K>
auto TensorExpression1<H, D, I>::operator*(TensorExpression3<F, D, J, I, K> const& other) {
   auto const lambda = [this, other](int j, int k) {
      using type = std::remove_cvref_t<decltype((*this)(0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         sum += (*this)(i) * other(j, i, k);
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, K>(lambda);
}

template<typename H, int D, char I>
template<typename F, char J, char K>
auto TensorExpression1<H, D, I>::operator*(TensorExpression3<F, D, J, K, I> const& other) {
   auto const lambda = [this, other](int j, int k) {
      using type = std::remove_cvref_t<decltype((*this)(0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         sum += (*this)(i) * other(j, k, i);
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, K>(lambda);
}

template<typename H, int D, char I>
template<typename F, char J, char K, char L>
auto TensorExpression1<H, D, I>::operator*(TensorExpression4<F, D, I, J, K, L> const& other) {
   auto const lambda = [this, other](int j, int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         sum += (*this)(i) * other(i, j, k, l);
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, K, L>(lambda);
}

template<typename H, int D, char I>
template<typename F, char J, char K, char L>
auto TensorExpression1<H, D, I>::operator*(TensorExpression4<F, D, J, I, K, L> const& other) {
   auto const lambda = [this, other](int j, int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         sum += (*this)(i) * other(j, i, k, l);
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, K, L>(lambda);
}

template<typename H, int D, char I>
template<typename F, char J, char K, char L>
auto TensorExpression1<H, D, I>::operator*(TensorExpression4<F, D, J, K, I, L> const& other) {
   auto const lambda = [this, other](int j, int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         sum += (*this)(i) * other(j, k, i, l);
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, K, L>(lambda);
}

template<typename H, int D, char I>
template<typename F, char J, char K, char L>
auto TensorExpression1<H, D, I>::operator*(TensorExpression4<F, D, J, K, L, I> const& other) {
   auto const lambda = [this, other](int j, int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         sum += (*this)(i) * other(j, k, l, i);
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, K, L>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression1<F, D, K> const& other) {
   auto const lambda = [this, other](int i, int j, int k) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0))>;
      auto sum = type(0);
      sum += (*this)(i, j) * other(k);
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, J, K>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression1<F, D, I> const& other) {
   auto const lambda = [this, other](int j) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         sum += (*this)(i, j) * other(i);
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, J>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression1<F, D, J> const& other) {
   auto const lambda = [this, other](int i) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         sum += (*this)(i, j) * other(j);
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, I>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K, char L>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression2<F, D, K, L> const& other) {
   auto const lambda = [this, other](int i, int j, int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0))>;
      auto sum = type(0);
      sum += (*this)(i, j) * other(k, l);
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, J, K, L>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression2<F, D, I, K> const& other) {
   auto const lambda = [this, other](int j, int k) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         sum += (*this)(i, j) * other(i, k);
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, K>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression2<F, D, K, I> const& other) {
   auto const lambda = [this, other](int j, int k) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         sum += (*this)(i, j) * other(k, i);
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, K>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression2<F, D, J, K> const& other) {
   auto const lambda = [this, other](int i, int k) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         sum += (*this)(i, j) * other(j, k);
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, K>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression2<F, D, K, J> const& other) {
   auto const lambda = [this, other](int i, int k) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         sum += (*this)(i, j) * other(k, j);
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, K>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression2<F, D, I, J> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j) * other(i, j);
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression2<F, D, J, I> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j) * other(j, i);
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K, char L>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression3<F, D, I, K, L> const& other) {
   auto const lambda = [this, other](int j, int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         sum += (*this)(i, j) * other(i, k, l);
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, K, L>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K, char L>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression3<F, D, K, I, L> const& other) {
   auto const lambda = [this, other](int j, int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         sum += (*this)(i, j) * other(k, i, l);
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, K, L>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K, char L>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression3<F, D, K, L, I> const& other) {
   auto const lambda = [this, other](int j, int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         sum += (*this)(i, j) * other(k, l, i);
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, K, L>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K, char L>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression3<F, D, J, K, L> const& other) {
   auto const lambda = [this, other](int i, int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         sum += (*this)(i, j) * other(j, k, l);
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, K, L>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K, char L>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression3<F, D, K, J, L> const& other) {
   auto const lambda = [this, other](int i, int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         sum += (*this)(i, j) * other(k, j, l);
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, K, L>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K, char L>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression3<F, D, K, L, J> const& other) {
   auto const lambda = [this, other](int i, int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         sum += (*this)(i, j) * other(k, l, j);
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, K, L>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression3<F, D, I, J, K> const& other) {
   auto const lambda = [this, other](int k) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j) * other(i, j, k);
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, K>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression3<F, D, I, K, J> const& other) {
   auto const lambda = [this, other](int k) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j) * other(i, k, j);
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, K>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression3<F, D, J, I, K> const& other) {
   auto const lambda = [this, other](int k) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j) * other(j, i, k);
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, K>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression3<F, D, J, K, I> const& other) {
   auto const lambda = [this, other](int k) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j) * other(j, k, i);
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, K>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression3<F, D, K, I, J> const& other) {
   auto const lambda = [this, other](int k) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j) * other(k, i, j);
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, K>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression3<F, D, K, J, I> const& other) {
   auto const lambda = [this, other](int k) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j) * other(k, j, i);
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, K>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K, char L, char M>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression4<F, D, I, K, L, M> const& other) {
   auto const lambda = [this, other](int j, int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         sum += (*this)(i, j) * other(i, k, l, m);
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, K, L, M>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K, char L, char M>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression4<F, D, K, I, L, M> const& other) {
   auto const lambda = [this, other](int j, int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         sum += (*this)(i, j) * other(k, i, l, m);
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, K, L, M>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K, char L, char M>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression4<F, D, K, L, I, M> const& other) {
   auto const lambda = [this, other](int j, int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         sum += (*this)(i, j) * other(k, l, i, m);
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, K, L, M>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K, char L, char M>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression4<F, D, K, L, M, I> const& other) {
   auto const lambda = [this, other](int j, int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         sum += (*this)(i, j) * other(k, l, m, i);
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, K, L, M>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K, char L, char M>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression4<F, D, J, K, L, M> const& other) {
   auto const lambda = [this, other](int i, int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         sum += (*this)(i, j) * other(j, k, l, m);
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, K, L, M>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K, char L, char M>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression4<F, D, K, J, L, M> const& other) {
   auto const lambda = [this, other](int i, int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         sum += (*this)(i, j) * other(k, j, l, m);
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, K, L, M>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K, char L, char M>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression4<F, D, K, L, J, M> const& other) {
   auto const lambda = [this, other](int i, int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         sum += (*this)(i, j) * other(k, l, j, m);
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, K, L, M>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K, char L, char M>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression4<F, D, K, L, M, J> const& other) {
   auto const lambda = [this, other](int i, int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         sum += (*this)(i, j) * other(k, l, m, j);
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, K, L, M>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K, char L>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression4<F, D, I, J, K, L> const& other) {
   auto const lambda = [this, other](int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j) * other(i, j, k, l);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, L>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K, char L>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression4<F, D, I, K, J, L> const& other) {
   auto const lambda = [this, other](int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j) * other(i, k, j, l);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, L>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K, char L>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression4<F, D, I, K, L, J> const& other) {
   auto const lambda = [this, other](int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j) * other(i, k, l, j);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, L>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K, char L>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression4<F, D, J, I, K, L> const& other) {
   auto const lambda = [this, other](int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j) * other(j, i, k, l);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, L>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K, char L>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression4<F, D, J, K, I, L> const& other) {
   auto const lambda = [this, other](int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j) * other(j, k, i, l);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, L>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K, char L>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression4<F, D, J, K, L, I> const& other) {
   auto const lambda = [this, other](int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j) * other(j, k, l, i);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, L>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K, char L>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression4<F, D, K, I, J, L> const& other) {
   auto const lambda = [this, other](int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j) * other(k, i, j, l);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, L>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K, char L>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression4<F, D, K, I, L, J> const& other) {
   auto const lambda = [this, other](int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j) * other(k, i, l, j);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, L>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K, char L>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression4<F, D, K, J, I, L> const& other) {
   auto const lambda = [this, other](int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j) * other(k, j, i, l);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, L>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K, char L>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression4<F, D, K, J, L, I> const& other) {
   auto const lambda = [this, other](int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j) * other(k, j, l, i);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, L>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K, char L>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression4<F, D, K, L, I, J> const& other) {
   auto const lambda = [this, other](int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j) * other(k, l, i, j);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, L>(lambda);
}

template<typename H, int D, char I, char J>
template<typename F, char K, char L>
auto TensorExpression2<H, D, I, J>::operator*(TensorExpression4<F, D, K, L, J, I> const& other) {
   auto const lambda = [this, other](int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j) * other(k, l, j, i);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression1<F, D, L> const& other) {
   auto const lambda = [this, other](int i, int j, int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0))>;
      auto sum = type(0);
      sum += (*this)(i, j, k) * other(l);
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, J, K, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression1<F, D, I> const& other) {
   auto const lambda = [this, other](int j, int k) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         sum += (*this)(i, j, k) * other(i);
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, K>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression1<F, D, J> const& other) {
   auto const lambda = [this, other](int i, int k) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         sum += (*this)(i, j, k) * other(j);
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, K>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression1<F, D, K> const& other) {
   auto const lambda = [this, other](int i, int j) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0))>;
      auto sum = type(0);
      for(int k = 0; k < D; k++ ) {
         sum += (*this)(i, j, k) * other(k);
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, J>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression2<F, D, I, L> const& other) {
   auto const lambda = [this, other](int j, int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         sum += (*this)(i, j, k) * other(i, l);
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, K, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression2<F, D, L, I> const& other) {
   auto const lambda = [this, other](int j, int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         sum += (*this)(i, j, k) * other(l, i);
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, K, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression2<F, D, J, L> const& other) {
   auto const lambda = [this, other](int i, int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         sum += (*this)(i, j, k) * other(j, l);
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, K, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression2<F, D, L, J> const& other) {
   auto const lambda = [this, other](int i, int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         sum += (*this)(i, j, k) * other(l, j);
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, K, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression2<F, D, I, J> const& other) {
   auto const lambda = [this, other](int k) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k) * other(i, j);
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, K>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression2<F, D, J, I> const& other) {
   auto const lambda = [this, other](int k) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k) * other(j, i);
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, K>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression2<F, D, K, L> const& other) {
   auto const lambda = [this, other](int i, int j, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int k = 0; k < D; k++ ) {
         sum += (*this)(i, j, k) * other(k, l);
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, J, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression2<F, D, L, K> const& other) {
   auto const lambda = [this, other](int i, int j, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int k = 0; k < D; k++ ) {
         sum += (*this)(i, j, k) * other(l, k);
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, J, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression2<F, D, I, K> const& other) {
   auto const lambda = [this, other](int j) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(i, k);
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, J>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression2<F, D, K, I> const& other) {
   auto const lambda = [this, other](int j) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(k, i);
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, J>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression2<F, D, J, K> const& other) {
   auto const lambda = [this, other](int i) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(j, k);
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, I>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression2<F, D, K, J> const& other) {
   auto const lambda = [this, other](int i) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(k, j);
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, I>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, I, L, M> const& other) {
   auto const lambda = [this, other](int j, int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         sum += (*this)(i, j, k) * other(i, l, m);
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, K, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, L, I, M> const& other) {
   auto const lambda = [this, other](int j, int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         sum += (*this)(i, j, k) * other(l, i, m);
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, K, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, L, M, I> const& other) {
   auto const lambda = [this, other](int j, int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         sum += (*this)(i, j, k) * other(l, m, i);
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, K, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, J, L, M> const& other) {
   auto const lambda = [this, other](int i, int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         sum += (*this)(i, j, k) * other(j, l, m);
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, K, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, L, J, M> const& other) {
   auto const lambda = [this, other](int i, int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         sum += (*this)(i, j, k) * other(l, j, m);
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, K, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, L, M, J> const& other) {
   auto const lambda = [this, other](int i, int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         sum += (*this)(i, j, k) * other(l, m, j);
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, K, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, I, J, L> const& other) {
   auto const lambda = [this, other](int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k) * other(i, j, l);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, I, L, J> const& other) {
   auto const lambda = [this, other](int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k) * other(i, l, j);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, J, I, L> const& other) {
   auto const lambda = [this, other](int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k) * other(j, i, l);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, J, L, I> const& other) {
   auto const lambda = [this, other](int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k) * other(j, l, i);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, L, I, J> const& other) {
   auto const lambda = [this, other](int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k) * other(l, i, j);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, L, J, I> const& other) {
   auto const lambda = [this, other](int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k) * other(l, j, i);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, K, L, M> const& other) {
   auto const lambda = [this, other](int i, int j, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int k = 0; k < D; k++ ) {
         sum += (*this)(i, j, k) * other(k, l, m);
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, J, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, L, K, M> const& other) {
   auto const lambda = [this, other](int i, int j, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int k = 0; k < D; k++ ) {
         sum += (*this)(i, j, k) * other(l, k, m);
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, J, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, L, M, K> const& other) {
   auto const lambda = [this, other](int i, int j, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int k = 0; k < D; k++ ) {
         sum += (*this)(i, j, k) * other(l, m, k);
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, J, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, I, K, L> const& other) {
   auto const lambda = [this, other](int j, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(i, k, l);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, I, L, K> const& other) {
   auto const lambda = [this, other](int j, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(i, l, k);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, K, I, L> const& other) {
   auto const lambda = [this, other](int j, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(k, i, l);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, K, L, I> const& other) {
   auto const lambda = [this, other](int j, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(k, l, i);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, L, I, K> const& other) {
   auto const lambda = [this, other](int j, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(l, i, k);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, L, K, I> const& other) {
   auto const lambda = [this, other](int j, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(l, k, i);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, J, K, L> const& other) {
   auto const lambda = [this, other](int i, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(j, k, l);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, J, L, K> const& other) {
   auto const lambda = [this, other](int i, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(j, l, k);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, K, J, L> const& other) {
   auto const lambda = [this, other](int i, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(k, j, l);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, K, L, J> const& other) {
   auto const lambda = [this, other](int i, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(k, l, j);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, L, J, K> const& other) {
   auto const lambda = [this, other](int i, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(l, j, k);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, L, K, J> const& other) {
   auto const lambda = [this, other](int i, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(l, k, j);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, I, J, K> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(i, j, k);
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, I, K, J> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(i, k, j);
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, J, I, K> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(j, i, k);
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, J, K, I> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(j, k, i);
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, K, I, J> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(k, i, j);
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression3<F, D, K, J, I> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(k, j, i);
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, I, J, L, M> const& other) {
   auto const lambda = [this, other](int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k) * other(i, j, l, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, K, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, I, L, J, M> const& other) {
   auto const lambda = [this, other](int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k) * other(i, l, j, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, K, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, I, L, M, J> const& other) {
   auto const lambda = [this, other](int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k) * other(i, l, m, j);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, K, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, J, I, L, M> const& other) {
   auto const lambda = [this, other](int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k) * other(j, i, l, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, K, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, J, L, I, M> const& other) {
   auto const lambda = [this, other](int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k) * other(j, l, i, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, K, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, J, L, M, I> const& other) {
   auto const lambda = [this, other](int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k) * other(j, l, m, i);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, K, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, L, I, J, M> const& other) {
   auto const lambda = [this, other](int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k) * other(l, i, j, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, K, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, L, I, M, J> const& other) {
   auto const lambda = [this, other](int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k) * other(l, i, m, j);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, K, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, L, J, I, M> const& other) {
   auto const lambda = [this, other](int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k) * other(l, j, i, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, K, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, L, J, M, I> const& other) {
   auto const lambda = [this, other](int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k) * other(l, j, m, i);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, K, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, L, M, I, J> const& other) {
   auto const lambda = [this, other](int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k) * other(l, m, i, j);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, K, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, L, M, J, I> const& other) {
   auto const lambda = [this, other](int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k) * other(l, m, j, i);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, K, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, I, K, L, M> const& other) {
   auto const lambda = [this, other](int j, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(i, k, l, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, I, L, K, M> const& other) {
   auto const lambda = [this, other](int j, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(i, l, k, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, I, L, M, K> const& other) {
   auto const lambda = [this, other](int j, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(i, l, m, k);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, K, I, L, M> const& other) {
   auto const lambda = [this, other](int j, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(k, i, l, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, K, L, I, M> const& other) {
   auto const lambda = [this, other](int j, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(k, l, i, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, K, L, M, I> const& other) {
   auto const lambda = [this, other](int j, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(k, l, m, i);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, L, I, K, M> const& other) {
   auto const lambda = [this, other](int j, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(l, i, k, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, L, I, M, K> const& other) {
   auto const lambda = [this, other](int j, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(l, i, m, k);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, L, K, I, M> const& other) {
   auto const lambda = [this, other](int j, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(l, k, i, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, L, K, M, I> const& other) {
   auto const lambda = [this, other](int j, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(l, k, m, i);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, L, M, I, K> const& other) {
   auto const lambda = [this, other](int j, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(l, m, i, k);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, L, M, K, I> const& other) {
   auto const lambda = [this, other](int j, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(l, m, k, i);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, J, K, L, M> const& other) {
   auto const lambda = [this, other](int i, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(j, k, l, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, J, L, K, M> const& other) {
   auto const lambda = [this, other](int i, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(j, l, k, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, J, L, M, K> const& other) {
   auto const lambda = [this, other](int i, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(j, l, m, k);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, K, J, L, M> const& other) {
   auto const lambda = [this, other](int i, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(k, j, l, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, K, L, J, M> const& other) {
   auto const lambda = [this, other](int i, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(k, l, j, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, K, L, M, J> const& other) {
   auto const lambda = [this, other](int i, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(k, l, m, j);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, L, J, K, M> const& other) {
   auto const lambda = [this, other](int i, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(l, j, k, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, L, J, M, K> const& other) {
   auto const lambda = [this, other](int i, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(l, j, m, k);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, L, K, J, M> const& other) {
   auto const lambda = [this, other](int i, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(l, k, j, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, L, K, M, J> const& other) {
   auto const lambda = [this, other](int i, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(l, k, m, j);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, L, M, J, K> const& other) {
   auto const lambda = [this, other](int i, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(l, m, j, k);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L, char M>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, L, M, K, J> const& other) {
   auto const lambda = [this, other](int i, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k) * other(l, m, k, j);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, I, J, K, L> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(i, j, k, l);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, I, J, L, K> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(i, j, l, k);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, I, K, J, L> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(i, k, j, l);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, I, K, L, J> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(i, k, l, j);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, I, L, J, K> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(i, l, j, k);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, I, L, K, J> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(i, l, k, j);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, J, I, K, L> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(j, i, k, l);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, J, I, L, K> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(j, i, l, k);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, J, K, I, L> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(j, k, i, l);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, J, K, L, I> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(j, k, l, i);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, J, L, I, K> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(j, l, i, k);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, J, L, K, I> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(j, l, k, i);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, K, I, J, L> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(k, i, j, l);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, K, I, L, J> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(k, i, l, j);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, K, J, I, L> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(k, j, i, l);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, K, J, L, I> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(k, j, l, i);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, K, L, I, J> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(k, l, i, j);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, K, L, J, I> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(k, l, j, i);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, L, I, J, K> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(l, i, j, k);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, L, I, K, J> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(l, i, k, j);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, L, J, I, K> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(l, j, i, k);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, L, J, K, I> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(l, j, k, i);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, L, K, I, J> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(l, k, i, j);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K>
template<typename F, char L>
auto TensorExpression3<H, D, I, J, K>::operator*(TensorExpression4<F, D, L, K, J, I> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k) * other(l, k, j, i);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression1<F, D, I> const& other) {
   auto const lambda = [this, other](int j, int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         sum += (*this)(i, j, k, l) * other(i);
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, K, L>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression1<F, D, J> const& other) {
   auto const lambda = [this, other](int i, int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         sum += (*this)(i, j, k, l) * other(j);
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, K, L>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression1<F, D, K> const& other) {
   auto const lambda = [this, other](int i, int j, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0))>;
      auto sum = type(0);
      for(int k = 0; k < D; k++ ) {
         sum += (*this)(i, j, k, l) * other(k);
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, J, L>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression1<F, D, L> const& other) {
   auto const lambda = [this, other](int i, int j, int k) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0))>;
      auto sum = type(0);
      for(int l = 0; l < D; l++ ) {
         sum += (*this)(i, j, k, l) * other(l);
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, J, K>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression2<F, D, I, M> const& other) {
   auto const lambda = [this, other](int j, int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         sum += (*this)(i, j, k, l) * other(i, m);
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, K, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression2<F, D, M, I> const& other) {
   auto const lambda = [this, other](int j, int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         sum += (*this)(i, j, k, l) * other(m, i);
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, K, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression2<F, D, J, M> const& other) {
   auto const lambda = [this, other](int i, int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         sum += (*this)(i, j, k, l) * other(j, m);
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, K, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression2<F, D, M, J> const& other) {
   auto const lambda = [this, other](int i, int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         sum += (*this)(i, j, k, l) * other(m, j);
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, K, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression2<F, D, I, J> const& other) {
   auto const lambda = [this, other](int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k, l) * other(i, j);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, L>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression2<F, D, J, I> const& other) {
   auto const lambda = [this, other](int k, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k, l) * other(j, i);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, L>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression2<F, D, K, M> const& other) {
   auto const lambda = [this, other](int i, int j, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int k = 0; k < D; k++ ) {
         sum += (*this)(i, j, k, l) * other(k, m);
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, J, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression2<F, D, M, K> const& other) {
   auto const lambda = [this, other](int i, int j, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int k = 0; k < D; k++ ) {
         sum += (*this)(i, j, k, l) * other(m, k);
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, J, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression2<F, D, I, K> const& other) {
   auto const lambda = [this, other](int j, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(i, k);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, L>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression2<F, D, K, I> const& other) {
   auto const lambda = [this, other](int j, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(k, i);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, L>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression2<F, D, J, K> const& other) {
   auto const lambda = [this, other](int i, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(j, k);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, L>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression2<F, D, K, J> const& other) {
   auto const lambda = [this, other](int i, int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(k, j);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, L>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression2<F, D, L, M> const& other) {
   auto const lambda = [this, other](int i, int j, int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int l = 0; l < D; l++ ) {
         sum += (*this)(i, j, k, l) * other(l, m);
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, J, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression2<F, D, M, L> const& other) {
   auto const lambda = [this, other](int i, int j, int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int l = 0; l < D; l++ ) {
         sum += (*this)(i, j, k, l) * other(m, l);
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, J, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression2<F, D, I, L> const& other) {
   auto const lambda = [this, other](int j, int k) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(i, l);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, K>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression2<F, D, L, I> const& other) {
   auto const lambda = [this, other](int j, int k) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(l, i);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, K>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression2<F, D, J, L> const& other) {
   auto const lambda = [this, other](int i, int k) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(j, l);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, K>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression2<F, D, L, J> const& other) {
   auto const lambda = [this, other](int i, int k) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(l, j);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, K>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression2<F, D, K, L> const& other) {
   auto const lambda = [this, other](int i, int j) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int k = 0; k < D; k++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(k, l);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, J>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression2<F, D, L, K> const& other) {
   auto const lambda = [this, other](int i, int j) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0))>;
      auto sum = type(0);
      for(int k = 0; k < D; k++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(l, k);
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, J>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, I, J, M> const& other) {
   auto const lambda = [this, other](int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k, l) * other(i, j, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, K, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, I, M, J> const& other) {
   auto const lambda = [this, other](int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k, l) * other(i, m, j);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, K, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, J, I, M> const& other) {
   auto const lambda = [this, other](int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k, l) * other(j, i, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, K, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, J, M, I> const& other) {
   auto const lambda = [this, other](int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k, l) * other(j, m, i);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, K, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, M, I, J> const& other) {
   auto const lambda = [this, other](int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k, l) * other(m, i, j);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, K, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, M, J, I> const& other) {
   auto const lambda = [this, other](int k, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k, l) * other(m, j, i);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, K, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, I, K, M> const& other) {
   auto const lambda = [this, other](int j, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(i, k, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, I, M, K> const& other) {
   auto const lambda = [this, other](int j, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(i, m, k);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, K, I, M> const& other) {
   auto const lambda = [this, other](int j, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(k, i, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, K, M, I> const& other) {
   auto const lambda = [this, other](int j, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(k, m, i);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, M, I, K> const& other) {
   auto const lambda = [this, other](int j, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(m, i, k);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, M, K, I> const& other) {
   auto const lambda = [this, other](int j, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(m, k, i);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, J, K, M> const& other) {
   auto const lambda = [this, other](int i, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(j, k, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, J, M, K> const& other) {
   auto const lambda = [this, other](int i, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(j, m, k);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, K, J, M> const& other) {
   auto const lambda = [this, other](int i, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(k, j, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, K, M, J> const& other) {
   auto const lambda = [this, other](int i, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(k, m, j);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, M, J, K> const& other) {
   auto const lambda = [this, other](int i, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(m, j, k);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, M, K, J> const& other) {
   auto const lambda = [this, other](int i, int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(m, k, j);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, I, J, K> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(i, j, k);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, I, K, J> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(i, k, j);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, J, I, K> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(j, i, k);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, J, K, I> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(j, k, i);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, K, I, J> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(k, i, j);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, K, J, I> const& other) {
   auto const lambda = [this, other](int l) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(k, j, i);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, L>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, I, L, M> const& other) {
   auto const lambda = [this, other](int j, int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(i, l, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, I, M, L> const& other) {
   auto const lambda = [this, other](int j, int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(i, m, l);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, L, I, M> const& other) {
   auto const lambda = [this, other](int j, int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(l, i, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, L, M, I> const& other) {
   auto const lambda = [this, other](int j, int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(l, m, i);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, M, I, L> const& other) {
   auto const lambda = [this, other](int j, int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(m, i, l);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, M, L, I> const& other) {
   auto const lambda = [this, other](int j, int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(m, l, i);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, J, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, J, L, M> const& other) {
   auto const lambda = [this, other](int i, int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(j, l, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, J, M, L> const& other) {
   auto const lambda = [this, other](int i, int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(j, m, l);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, L, J, M> const& other) {
   auto const lambda = [this, other](int i, int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(l, j, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, L, M, J> const& other) {
   auto const lambda = [this, other](int i, int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(l, m, j);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, M, J, L> const& other) {
   auto const lambda = [this, other](int i, int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(m, j, l);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, M, L, J> const& other) {
   auto const lambda = [this, other](int i, int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(m, l, j);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, I, J, L> const& other) {
   auto const lambda = [this, other](int k) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(i, j, l);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, K>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, I, L, J> const& other) {
   auto const lambda = [this, other](int k) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(i, l, j);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, K>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, J, I, L> const& other) {
   auto const lambda = [this, other](int k) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(j, i, l);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, K>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, J, L, I> const& other) {
   auto const lambda = [this, other](int k) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(j, l, i);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, K>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, L, I, J> const& other) {
   auto const lambda = [this, other](int k) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(l, i, j);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, K>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, L, J, I> const& other) {
   auto const lambda = [this, other](int k) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(l, j, i);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, K>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, K, L, M> const& other) {
   auto const lambda = [this, other](int i, int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int k = 0; k < D; k++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(k, l, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, K, M, L> const& other) {
   auto const lambda = [this, other](int i, int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int k = 0; k < D; k++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(k, m, l);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, L, K, M> const& other) {
   auto const lambda = [this, other](int i, int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int k = 0; k < D; k++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(l, k, m);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, L, M, K> const& other) {
   auto const lambda = [this, other](int i, int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int k = 0; k < D; k++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(l, m, k);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, M, K, L> const& other) {
   auto const lambda = [this, other](int i, int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int k = 0; k < D; k++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(m, k, l);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, M, L, K> const& other) {
   auto const lambda = [this, other](int i, int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int k = 0; k < D; k++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(m, l, k);
         }
      }
      return sum;
   };
   return Tensor3Expression<decltype(lambda), D, I, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, I, K, L> const& other) {
   auto const lambda = [this, other](int j) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(i, k, l);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, J>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, I, L, K> const& other) {
   auto const lambda = [this, other](int j) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(i, l, k);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, J>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, K, I, L> const& other) {
   auto const lambda = [this, other](int j) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(k, i, l);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, J>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, K, L, I> const& other) {
   auto const lambda = [this, other](int j) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(k, l, i);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, J>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, L, I, K> const& other) {
   auto const lambda = [this, other](int j) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(l, i, k);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, J>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, L, K, I> const& other) {
   auto const lambda = [this, other](int j) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(l, k, i);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, J>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, J, K, L> const& other) {
   auto const lambda = [this, other](int i) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(j, k, l);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, I>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, J, L, K> const& other) {
   auto const lambda = [this, other](int i) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(j, l, k);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, I>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, K, J, L> const& other) {
   auto const lambda = [this, other](int i) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(k, j, l);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, I>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, K, L, J> const& other) {
   auto const lambda = [this, other](int i) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(k, l, j);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, I>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, L, J, K> const& other) {
   auto const lambda = [this, other](int i) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(l, j, k);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, I>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression3<F, D, L, K, J> const& other) {
   auto const lambda = [this, other](int i) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(l, k, j);
            }
         }
      }
      return sum;
   };
   return Tensor1Expression<decltype(lambda), D, I>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, J, M, N> const& other) {
   auto const lambda = [this, other](int k, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k, l) * other(i, j, m, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, K, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, M, J, N> const& other) {
   auto const lambda = [this, other](int k, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k, l) * other(i, m, j, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, K, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, M, N, J> const& other) {
   auto const lambda = [this, other](int k, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k, l) * other(i, m, n, j);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, K, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, I, M, N> const& other) {
   auto const lambda = [this, other](int k, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k, l) * other(j, i, m, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, K, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, M, I, N> const& other) {
   auto const lambda = [this, other](int k, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k, l) * other(j, m, i, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, K, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, M, N, I> const& other) {
   auto const lambda = [this, other](int k, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k, l) * other(j, m, n, i);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, K, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, I, J, N> const& other) {
   auto const lambda = [this, other](int k, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k, l) * other(m, i, j, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, K, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, I, N, J> const& other) {
   auto const lambda = [this, other](int k, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k, l) * other(m, i, n, j);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, K, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, J, I, N> const& other) {
   auto const lambda = [this, other](int k, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k, l) * other(m, j, i, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, K, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, J, N, I> const& other) {
   auto const lambda = [this, other](int k, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k, l) * other(m, j, n, i);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, K, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, N, I, J> const& other) {
   auto const lambda = [this, other](int k, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k, l) * other(m, n, i, j);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, K, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, N, J, I> const& other) {
   auto const lambda = [this, other](int k, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            sum += (*this)(i, j, k, l) * other(m, n, j, i);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, K, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, K, M, N> const& other) {
   auto const lambda = [this, other](int j, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(i, k, m, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, M, K, N> const& other) {
   auto const lambda = [this, other](int j, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(i, m, k, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, M, N, K> const& other) {
   auto const lambda = [this, other](int j, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(i, m, n, k);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, I, M, N> const& other) {
   auto const lambda = [this, other](int j, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(k, i, m, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, M, I, N> const& other) {
   auto const lambda = [this, other](int j, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(k, m, i, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, M, N, I> const& other) {
   auto const lambda = [this, other](int j, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(k, m, n, i);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, I, K, N> const& other) {
   auto const lambda = [this, other](int j, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(m, i, k, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, I, N, K> const& other) {
   auto const lambda = [this, other](int j, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(m, i, n, k);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, K, I, N> const& other) {
   auto const lambda = [this, other](int j, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(m, k, i, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, K, N, I> const& other) {
   auto const lambda = [this, other](int j, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(m, k, n, i);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, N, I, K> const& other) {
   auto const lambda = [this, other](int j, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(m, n, i, k);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, N, K, I> const& other) {
   auto const lambda = [this, other](int j, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(m, n, k, i);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, K, M, N> const& other) {
   auto const lambda = [this, other](int i, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(j, k, m, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, M, K, N> const& other) {
   auto const lambda = [this, other](int i, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(j, m, k, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, M, N, K> const& other) {
   auto const lambda = [this, other](int i, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(j, m, n, k);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, J, M, N> const& other) {
   auto const lambda = [this, other](int i, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(k, j, m, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, M, J, N> const& other) {
   auto const lambda = [this, other](int i, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(k, m, j, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, M, N, J> const& other) {
   auto const lambda = [this, other](int i, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(k, m, n, j);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, J, K, N> const& other) {
   auto const lambda = [this, other](int i, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(m, j, k, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, J, N, K> const& other) {
   auto const lambda = [this, other](int i, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(m, j, n, k);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, K, J, N> const& other) {
   auto const lambda = [this, other](int i, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(m, k, j, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, K, N, J> const& other) {
   auto const lambda = [this, other](int i, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(m, k, n, j);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, N, J, K> const& other) {
   auto const lambda = [this, other](int i, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(m, n, j, k);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, N, K, J> const& other) {
   auto const lambda = [this, other](int i, int l, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            sum += (*this)(i, j, k, l) * other(m, n, k, j);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, L, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, J, K, M> const& other) {
   auto const lambda = [this, other](int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(i, j, k, m);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, J, M, K> const& other) {
   auto const lambda = [this, other](int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(i, j, m, k);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, K, J, M> const& other) {
   auto const lambda = [this, other](int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(i, k, j, m);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, K, M, J> const& other) {
   auto const lambda = [this, other](int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(i, k, m, j);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, M, J, K> const& other) {
   auto const lambda = [this, other](int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(i, m, j, k);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, M, K, J> const& other) {
   auto const lambda = [this, other](int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(i, m, k, j);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, I, K, M> const& other) {
   auto const lambda = [this, other](int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(j, i, k, m);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, I, M, K> const& other) {
   auto const lambda = [this, other](int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(j, i, m, k);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, K, I, M> const& other) {
   auto const lambda = [this, other](int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(j, k, i, m);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, K, M, I> const& other) {
   auto const lambda = [this, other](int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(j, k, m, i);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, M, I, K> const& other) {
   auto const lambda = [this, other](int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(j, m, i, k);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, M, K, I> const& other) {
   auto const lambda = [this, other](int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(j, m, k, i);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, I, J, M> const& other) {
   auto const lambda = [this, other](int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(k, i, j, m);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, I, M, J> const& other) {
   auto const lambda = [this, other](int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(k, i, m, j);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, J, I, M> const& other) {
   auto const lambda = [this, other](int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(k, j, i, m);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, J, M, I> const& other) {
   auto const lambda = [this, other](int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(k, j, m, i);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, M, I, J> const& other) {
   auto const lambda = [this, other](int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(k, m, i, j);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, M, J, I> const& other) {
   auto const lambda = [this, other](int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(k, m, j, i);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, I, J, K> const& other) {
   auto const lambda = [this, other](int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(m, i, j, k);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, I, K, J> const& other) {
   auto const lambda = [this, other](int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(m, i, k, j);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, J, I, K> const& other) {
   auto const lambda = [this, other](int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(m, j, i, k);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, J, K, I> const& other) {
   auto const lambda = [this, other](int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(m, j, k, i);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, K, I, J> const& other) {
   auto const lambda = [this, other](int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(m, k, i, j);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, K, J, I> const& other) {
   auto const lambda = [this, other](int l, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               sum += (*this)(i, j, k, l) * other(m, k, j, i);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, L, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, L, M, N> const& other) {
   auto const lambda = [this, other](int j, int k, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(i, l, m, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, K, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, M, L, N> const& other) {
   auto const lambda = [this, other](int j, int k, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(i, m, l, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, K, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, M, N, L> const& other) {
   auto const lambda = [this, other](int j, int k, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(i, m, n, l);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, K, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, I, M, N> const& other) {
   auto const lambda = [this, other](int j, int k, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(l, i, m, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, K, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, M, I, N> const& other) {
   auto const lambda = [this, other](int j, int k, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(l, m, i, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, K, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, M, N, I> const& other) {
   auto const lambda = [this, other](int j, int k, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(l, m, n, i);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, K, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, I, L, N> const& other) {
   auto const lambda = [this, other](int j, int k, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(m, i, l, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, K, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, I, N, L> const& other) {
   auto const lambda = [this, other](int j, int k, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(m, i, n, l);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, K, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, L, I, N> const& other) {
   auto const lambda = [this, other](int j, int k, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(m, l, i, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, K, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, L, N, I> const& other) {
   auto const lambda = [this, other](int j, int k, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(m, l, n, i);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, K, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, N, I, L> const& other) {
   auto const lambda = [this, other](int j, int k, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(m, n, i, l);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, K, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, N, L, I> const& other) {
   auto const lambda = [this, other](int j, int k, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(m, n, l, i);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, J, K, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, L, M, N> const& other) {
   auto const lambda = [this, other](int i, int k, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(j, l, m, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, K, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, M, L, N> const& other) {
   auto const lambda = [this, other](int i, int k, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(j, m, l, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, K, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, M, N, L> const& other) {
   auto const lambda = [this, other](int i, int k, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(j, m, n, l);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, K, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, J, M, N> const& other) {
   auto const lambda = [this, other](int i, int k, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(l, j, m, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, K, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, M, J, N> const& other) {
   auto const lambda = [this, other](int i, int k, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(l, m, j, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, K, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, M, N, J> const& other) {
   auto const lambda = [this, other](int i, int k, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(l, m, n, j);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, K, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, J, L, N> const& other) {
   auto const lambda = [this, other](int i, int k, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(m, j, l, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, K, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, J, N, L> const& other) {
   auto const lambda = [this, other](int i, int k, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(m, j, n, l);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, K, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, L, J, N> const& other) {
   auto const lambda = [this, other](int i, int k, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(m, l, j, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, K, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, L, N, J> const& other) {
   auto const lambda = [this, other](int i, int k, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(m, l, n, j);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, K, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, N, J, L> const& other) {
   auto const lambda = [this, other](int i, int k, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(m, n, j, l);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, K, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, N, L, J> const& other) {
   auto const lambda = [this, other](int i, int k, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(m, n, l, j);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, K, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, J, L, M> const& other) {
   auto const lambda = [this, other](int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(i, j, l, m);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, J, M, L> const& other) {
   auto const lambda = [this, other](int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(i, j, m, l);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, L, J, M> const& other) {
   auto const lambda = [this, other](int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(i, l, j, m);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, L, M, J> const& other) {
   auto const lambda = [this, other](int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(i, l, m, j);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, M, J, L> const& other) {
   auto const lambda = [this, other](int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(i, m, j, l);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, M, L, J> const& other) {
   auto const lambda = [this, other](int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(i, m, l, j);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, I, L, M> const& other) {
   auto const lambda = [this, other](int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(j, i, l, m);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, I, M, L> const& other) {
   auto const lambda = [this, other](int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(j, i, m, l);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, L, I, M> const& other) {
   auto const lambda = [this, other](int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(j, l, i, m);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, L, M, I> const& other) {
   auto const lambda = [this, other](int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(j, l, m, i);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, M, I, L> const& other) {
   auto const lambda = [this, other](int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(j, m, i, l);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, M, L, I> const& other) {
   auto const lambda = [this, other](int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(j, m, l, i);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, I, J, M> const& other) {
   auto const lambda = [this, other](int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(l, i, j, m);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, I, M, J> const& other) {
   auto const lambda = [this, other](int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(l, i, m, j);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, J, I, M> const& other) {
   auto const lambda = [this, other](int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(l, j, i, m);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, J, M, I> const& other) {
   auto const lambda = [this, other](int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(l, j, m, i);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, M, I, J> const& other) {
   auto const lambda = [this, other](int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(l, m, i, j);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, M, J, I> const& other) {
   auto const lambda = [this, other](int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(l, m, j, i);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, I, J, L> const& other) {
   auto const lambda = [this, other](int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(m, i, j, l);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, I, L, J> const& other) {
   auto const lambda = [this, other](int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(m, i, l, j);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, J, I, L> const& other) {
   auto const lambda = [this, other](int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(m, j, i, l);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, J, L, I> const& other) {
   auto const lambda = [this, other](int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(m, j, l, i);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, L, I, J> const& other) {
   auto const lambda = [this, other](int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(m, l, i, j);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, L, J, I> const& other) {
   auto const lambda = [this, other](int k, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(m, l, j, i);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, K, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, L, M, N> const& other) {
   auto const lambda = [this, other](int i, int j, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int k = 0; k < D; k++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(k, l, m, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, J, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, M, L, N> const& other) {
   auto const lambda = [this, other](int i, int j, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int k = 0; k < D; k++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(k, m, l, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, J, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, M, N, L> const& other) {
   auto const lambda = [this, other](int i, int j, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int k = 0; k < D; k++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(k, m, n, l);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, J, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, K, M, N> const& other) {
   auto const lambda = [this, other](int i, int j, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int k = 0; k < D; k++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(l, k, m, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, J, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, M, K, N> const& other) {
   auto const lambda = [this, other](int i, int j, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int k = 0; k < D; k++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(l, m, k, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, J, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, M, N, K> const& other) {
   auto const lambda = [this, other](int i, int j, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int k = 0; k < D; k++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(l, m, n, k);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, J, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, K, L, N> const& other) {
   auto const lambda = [this, other](int i, int j, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int k = 0; k < D; k++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(m, k, l, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, J, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, K, N, L> const& other) {
   auto const lambda = [this, other](int i, int j, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int k = 0; k < D; k++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(m, k, n, l);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, J, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, L, K, N> const& other) {
   auto const lambda = [this, other](int i, int j, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int k = 0; k < D; k++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(m, l, k, n);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, J, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, L, N, K> const& other) {
   auto const lambda = [this, other](int i, int j, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int k = 0; k < D; k++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(m, l, n, k);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, J, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, N, K, L> const& other) {
   auto const lambda = [this, other](int i, int j, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int k = 0; k < D; k++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(m, n, k, l);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, J, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M, char N>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, N, L, K> const& other) {
   auto const lambda = [this, other](int i, int j, int m, int n) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int k = 0; k < D; k++ ) {
         for(int l = 0; l < D; l++ ) {
            sum += (*this)(i, j, k, l) * other(m, n, l, k);
         }
      }
      return sum;
   };
   return Tensor4Expression<decltype(lambda), D, I, J, M, N>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, K, L, M> const& other) {
   auto const lambda = [this, other](int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(i, k, l, m);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, K, M, L> const& other) {
   auto const lambda = [this, other](int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(i, k, m, l);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, L, K, M> const& other) {
   auto const lambda = [this, other](int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(i, l, k, m);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, L, M, K> const& other) {
   auto const lambda = [this, other](int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(i, l, m, k);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, M, K, L> const& other) {
   auto const lambda = [this, other](int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(i, m, k, l);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, M, L, K> const& other) {
   auto const lambda = [this, other](int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(i, m, l, k);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, I, L, M> const& other) {
   auto const lambda = [this, other](int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(k, i, l, m);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, I, M, L> const& other) {
   auto const lambda = [this, other](int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(k, i, m, l);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, L, I, M> const& other) {
   auto const lambda = [this, other](int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(k, l, i, m);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, L, M, I> const& other) {
   auto const lambda = [this, other](int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(k, l, m, i);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, M, I, L> const& other) {
   auto const lambda = [this, other](int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(k, m, i, l);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, M, L, I> const& other) {
   auto const lambda = [this, other](int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(k, m, l, i);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, I, K, M> const& other) {
   auto const lambda = [this, other](int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(l, i, k, m);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, I, M, K> const& other) {
   auto const lambda = [this, other](int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(l, i, m, k);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, K, I, M> const& other) {
   auto const lambda = [this, other](int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(l, k, i, m);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, K, M, I> const& other) {
   auto const lambda = [this, other](int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(l, k, m, i);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, M, I, K> const& other) {
   auto const lambda = [this, other](int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(l, m, i, k);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, M, K, I> const& other) {
   auto const lambda = [this, other](int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(l, m, k, i);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, I, K, L> const& other) {
   auto const lambda = [this, other](int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(m, i, k, l);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, I, L, K> const& other) {
   auto const lambda = [this, other](int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(m, i, l, k);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, K, I, L> const& other) {
   auto const lambda = [this, other](int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(m, k, i, l);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, K, L, I> const& other) {
   auto const lambda = [this, other](int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(m, k, l, i);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, L, I, K> const& other) {
   auto const lambda = [this, other](int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(m, l, i, k);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, L, K, I> const& other) {
   auto const lambda = [this, other](int j, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(m, l, k, i);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, J, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, K, L, M> const& other) {
   auto const lambda = [this, other](int i, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(j, k, l, m);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, K, M, L> const& other) {
   auto const lambda = [this, other](int i, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(j, k, m, l);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, L, K, M> const& other) {
   auto const lambda = [this, other](int i, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(j, l, k, m);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, L, M, K> const& other) {
   auto const lambda = [this, other](int i, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(j, l, m, k);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, M, K, L> const& other) {
   auto const lambda = [this, other](int i, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(j, m, k, l);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, M, L, K> const& other) {
   auto const lambda = [this, other](int i, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(j, m, l, k);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, J, L, M> const& other) {
   auto const lambda = [this, other](int i, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(k, j, l, m);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, J, M, L> const& other) {
   auto const lambda = [this, other](int i, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(k, j, m, l);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, L, J, M> const& other) {
   auto const lambda = [this, other](int i, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(k, l, j, m);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, L, M, J> const& other) {
   auto const lambda = [this, other](int i, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(k, l, m, j);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, M, J, L> const& other) {
   auto const lambda = [this, other](int i, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(k, m, j, l);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, M, L, J> const& other) {
   auto const lambda = [this, other](int i, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(k, m, l, j);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, J, K, M> const& other) {
   auto const lambda = [this, other](int i, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(l, j, k, m);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, J, M, K> const& other) {
   auto const lambda = [this, other](int i, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(l, j, m, k);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, K, J, M> const& other) {
   auto const lambda = [this, other](int i, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(l, k, j, m);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, K, M, J> const& other) {
   auto const lambda = [this, other](int i, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(l, k, m, j);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, M, J, K> const& other) {
   auto const lambda = [this, other](int i, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(l, m, j, k);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, M, K, J> const& other) {
   auto const lambda = [this, other](int i, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(l, m, k, j);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, J, K, L> const& other) {
   auto const lambda = [this, other](int i, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(m, j, k, l);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, J, L, K> const& other) {
   auto const lambda = [this, other](int i, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(m, j, l, k);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, K, J, L> const& other) {
   auto const lambda = [this, other](int i, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(m, k, j, l);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, K, L, J> const& other) {
   auto const lambda = [this, other](int i, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(m, k, l, j);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, L, J, K> const& other) {
   auto const lambda = [this, other](int i, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(m, l, j, k);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F, char M>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, M, L, K, J> const& other) {
   auto const lambda = [this, other](int i, int m) {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int j = 0; j < D; j++ ) {
         for(int k = 0; k < D; k++ ) {
            for(int l = 0; l < D; l++ ) {
               sum += (*this)(i, j, k, l) * other(m, l, k, j);
            }
         }
      }
      return sum;
   };
   return Tensor2Expression<decltype(lambda), D, I, M>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, J, K, L> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               for(int l = 0; l < D; l++ ) {
                  sum += (*this)(i, j, k, l) * other(i, j, k, l);
               }
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, J, L, K> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               for(int l = 0; l < D; l++ ) {
                  sum += (*this)(i, j, k, l) * other(i, j, l, k);
               }
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, K, J, L> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               for(int l = 0; l < D; l++ ) {
                  sum += (*this)(i, j, k, l) * other(i, k, j, l);
               }
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, K, L, J> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               for(int l = 0; l < D; l++ ) {
                  sum += (*this)(i, j, k, l) * other(i, k, l, j);
               }
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, L, J, K> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               for(int l = 0; l < D; l++ ) {
                  sum += (*this)(i, j, k, l) * other(i, l, j, k);
               }
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, I, L, K, J> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               for(int l = 0; l < D; l++ ) {
                  sum += (*this)(i, j, k, l) * other(i, l, k, j);
               }
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, I, K, L> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               for(int l = 0; l < D; l++ ) {
                  sum += (*this)(i, j, k, l) * other(j, i, k, l);
               }
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, I, L, K> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               for(int l = 0; l < D; l++ ) {
                  sum += (*this)(i, j, k, l) * other(j, i, l, k);
               }
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, K, I, L> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               for(int l = 0; l < D; l++ ) {
                  sum += (*this)(i, j, k, l) * other(j, k, i, l);
               }
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, K, L, I> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               for(int l = 0; l < D; l++ ) {
                  sum += (*this)(i, j, k, l) * other(j, k, l, i);
               }
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, L, I, K> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               for(int l = 0; l < D; l++ ) {
                  sum += (*this)(i, j, k, l) * other(j, l, i, k);
               }
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, J, L, K, I> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               for(int l = 0; l < D; l++ ) {
                  sum += (*this)(i, j, k, l) * other(j, l, k, i);
               }
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, I, J, L> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               for(int l = 0; l < D; l++ ) {
                  sum += (*this)(i, j, k, l) * other(k, i, j, l);
               }
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, I, L, J> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               for(int l = 0; l < D; l++ ) {
                  sum += (*this)(i, j, k, l) * other(k, i, l, j);
               }
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, J, I, L> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               for(int l = 0; l < D; l++ ) {
                  sum += (*this)(i, j, k, l) * other(k, j, i, l);
               }
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, J, L, I> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               for(int l = 0; l < D; l++ ) {
                  sum += (*this)(i, j, k, l) * other(k, j, l, i);
               }
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, L, I, J> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               for(int l = 0; l < D; l++ ) {
                  sum += (*this)(i, j, k, l) * other(k, l, i, j);
               }
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, K, L, J, I> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               for(int l = 0; l < D; l++ ) {
                  sum += (*this)(i, j, k, l) * other(k, l, j, i);
               }
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, I, J, K> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               for(int l = 0; l < D; l++ ) {
                  sum += (*this)(i, j, k, l) * other(l, i, j, k);
               }
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, I, K, J> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               for(int l = 0; l < D; l++ ) {
                  sum += (*this)(i, j, k, l) * other(l, i, k, j);
               }
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, J, I, K> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               for(int l = 0; l < D; l++ ) {
                  sum += (*this)(i, j, k, l) * other(l, j, i, k);
               }
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, J, K, I> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               for(int l = 0; l < D; l++ ) {
                  sum += (*this)(i, j, k, l) * other(l, j, k, i);
               }
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, K, I, J> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               for(int l = 0; l < D; l++ ) {
                  sum += (*this)(i, j, k, l) * other(l, k, i, j);
               }
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

template<typename H, int D, char I, char J, char K, char L>
template<typename F>
auto TensorExpression4<H, D, I, J, K, L>::operator*(TensorExpression4<F, D, L, K, J, I> const& other) {
   auto const lambda = [this, other]() {
      using type = std::remove_cvref_t<decltype((*this)(0, 0, 0, 0) * other(0, 0, 0, 0))>;
      auto sum = type(0);
      for(int i = 0; i < D; i++ ) {
         for(int j = 0; j < D; j++ ) {
            for(int k = 0; k < D; k++ ) {
               for(int l = 0; l < D; l++ ) {
                  sum += (*this)(i, j, k, l) * other(l, k, j, i);
               }
            }
         }
      }
      return sum;
   };
   return Tensor0Expression<decltype(lambda), D>(lambda);
}

