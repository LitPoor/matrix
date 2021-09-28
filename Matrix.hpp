//
// Created by litpoor on 2021/9/26.
//

#ifndef EXERCISE_MATRIX_HPP
#define EXERCISE_MATRIX_HPP

#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;

template<typename T>
class Matrix
{

private:
    int rows, cols;
    vector<vector<T>> matrix;

public:

    Matrix(int rows, int cols, T e = 0); //构造一个元素为e的rows×cols矩阵

    Matrix(vector<vector<T>> &vv); //用vector<vector<T>> 构造一个矩阵

    void show() const;

    int det();  // 求矩阵行列式

    Matrix companion();      // 求伴随矩阵

    Matrix inverse();      // 求逆矩阵

    T getElement(int rows, int cols);   //获取矩阵中元素


    int getCols()
    { return cols; }

    int getRows()
    { return rows; }

    /// 如果把模板类的友元函数的定义放在类外会报错，加上template相当于重新定义了一个模板函数
    friend ostream &operator<<(ostream &o, Matrix<T> &m)
    {
        m.show();
        return o;
    }

    template<typename K>
    bool operator==(Matrix<K> &m) const;


    template<typename K>
    Matrix<K>& operator=(Matrix<K> &m);

    template<typename K, typename V>
    Matrix<V> operator+(Matrix<K> &m);


};


/*
 * 用类似归并排序的方法求逆序对的数量,时间复杂度o(nlogn)
 */
int merge_sort(vector<int> &v, int l, int r)
{
    if (l >= r) return 0;

    vector<int> temp;

    int mid = l + r >> 1;   //分解子问题
    int res = merge_sort(v, l, mid) + merge_sort(v, mid + 1, r);  //递归处理子问题
    int k = 0, i = l, j = mid + 1;

    //合并子问题
    while (i <= mid && j <= r)
    {
        if (v[i] <= v[j]) temp.push_back(v[i++]);
        else
        {
            temp.push_back(v[j++]);
            res += mid - i + 1;
        }
    }
    while (i <= mid) temp.push_back(v[i++]);
    while (j <= l) temp.push_back(v[j++]);
    for (int i = 0; i < temp.size(); i++)
        v[i] = temp[i];

    return res;
}


template<typename T>
Matrix<T>::Matrix(vector<vector<T>> &vv):matrix(vv)
{
    rows = vv.size();
    cols = vv[0].size();
}


template<typename T>
Matrix<T>::Matrix(int rows, int cols, T e):rows(rows), cols(cols)
{
    for (int i = 0; i < rows; i++)
    {
        vector<T> v;
        for (int j = 0; j < cols; j++)
        {
            v.push_back(e);
        }
        matrix.push_back(v);
    }
}


template<typename T>
void Matrix<T>::show() const
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
            cout << matrix[i][j] << ' ';
        cout << "\n";
    }
}


/*
 * 用行列式的定义求行列式
 * 设N阶矩阵每列元素的列下标分别为0，1，2..N，求出该序列的全排列，并用vector<vector<int>>储存，vector<int>中存储单个排序
 * 计算每个排列的逆序对数量，算出其所对应元素的前缀为1还是-1
 * 分别以0，1，2...N为元素所对应的行下标，vector<int>中所存的值为列下标，存成元素的下标，然后相乘
 * 最后求和
 * 时间复杂度o(n!)   っ(。>︿<)_θ
 */
template<typename T>
int Matrix<T>::det()
{
    if (cols != rows)
    {
        cout << "Error: 只能求方阵的行列式哦！\n";
        exit(0);
    }
    vector<int> v;
    for (int i = 0; i < rows; i++) v.push_back(i);   //生成一个顺序序列
    vector<vector<int>> vv; //用来存放该序列的全排列集合
    int ans = 0;   //计算结果

    do
    {
        vv.push_back(v);
    } while (next_permutation(v.begin(), v.end()));  //STL的一个模板函数，可以取得所标示之序列的下一个排列组合

    for (int i = 0; i < vv.size(); i++)
    {

        int res = 1;
        for (int j = 0; j < rows; j++)
            res *= matrix[j][vv[i][j]];           //计算元素乘积

        int p = merge_sort(vv[i], 0, rows - 1);
        if (p % 2) res *= -1;  //计算前缀的正负

        ans += res;
    }


    return ans;
}


template<typename T>
T Matrix<T>::getElement(int rows, int cols)
{
    if (rows > this->rows || cols > this->cols)
    {
        cout << "Error: 超出范围\n";
        exit(0);
    }
    return matrix[rows - 1][cols - 1];
}


template<typename T>
template<typename K>
bool Matrix<T>::operator==(Matrix<K> &m) const
{
    if (m.rows != this->rows || m.cols != this->cols)
        return false;
    else
    {
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                if (this->matrix[i][j] != m.matrix[i][j])
                    return false;
        return true;
    }
}


template<typename T>
template<typename K>
Matrix<K>& Matrix<T>::operator=(Matrix<K> &m)
{
    this->rows = m.rows;
    this->cols = m.cols;
    this->matrix = m.matrix;
    return *this;
}


template<typename T>
template<typename K, typename V>
Matrix<V> Matrix<T>::operator+(Matrix<K> &m)
{
    if (m.rows != this->rows || m.cols != this->cols)
    {
        cout << "Error: 用于矩阵加法的维度不对，请检查并确保第一个矩阵的行和列数于第二个矩阵相等\n";
        exit(0);
    }

    string s1 = typeid(m.matrix[0][0]).name();  //获取矩阵的数据类型，该返回值为一个字符串常量
    string s2 = typeid(this->matrix[0][0]).name();


    ///基本数据类型的返回值,测了半天qwq
    ///int -> i   short -> s    long -> l   long long -> x
    ///float -> f  long double -> e   double -> d

    if ((s1 != "i" && s1 != "s" && s1 != "l" && s1 != "x" && s1 != "f" && s1 != "e" && s1 != "d")
        || (s2 != "i" && s2 != "s" && s2 != "l" && s2 != "x" && s2 != "f" && s2 != "e" && s2 != "d"))
    {
        cout << "Error: 矩阵的运算只支持整形和浮点型\n";
        exit(0);
    }


    ///如果运算的矩阵有一个类型为浮点型，则返回矩阵的类型为double
    if (s1 == "f" || s1 == "e" || s1 == "d" || s2 == "f" || s2 == "e" || s2 == "d")
    {
        vector<double> v(cols);
        vector<vector<double>> vv;

        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
                v[j] = this->matrix[i][j] + m.matrix[i][j];

            vv.push_back(v);
        }

        return Matrix<double>(vv);
    } else if (s1 == "x" || s2 == "x")   //如果运算的矩阵有一个类型为long long，则返回矩阵的类型为long long
    {
        vector<long long> v(cols);
        vector<vector<long long>> vv;

        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
                v[j] = this->matrix[i][j] + m.matrix[i][j];

            vv.push_back(v);
        }

        return Matrix<long long>(vv);
    } else
    {
        vector<int> v(cols);
        vector<vector<int>> vv;

        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
                v[j] = this->matrix[i][j] + m.matrix[i][j];

            vv.push_back(v);
        }

        return Matrix<int>(vv);
    }

}


#endif //EXERCISE_MATRIX_HPP
