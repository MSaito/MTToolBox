;; -*- coding:utf-8 -*-

XORSHIFT のサンプル

このサンプルディレクトリは XORSHIFT[1](2003 Marsaglia) の
均等分布次元をMTToolBox を使用して求めるサンプルである。
同時に、MTToolBox 全体のチュートリアルとなっている。

1. xorshift 128ビット版の最小多項式を求める。
まず最初に周期が最大であることを確認する。
xorshift-1.cpp では、xorshift の最小多項式を求めて原始性を判定している。
シフト量などのパラメータは[1]より採ったもの。

2. xorshift 128 ビット版のvビット精度均等分布次元k(v)を求める。
xorshift-2.cpp では xorshift のvビット精度均等分布次元を求めている。

3. xorshift 128 ビット版のシフト量を変えてk(v)を求める。
ただし、最大周期を達成しなければk(v)は計算しない。
出力をソートすればdeltaを最小にするパラメータがわかる。
0から31までのシフト量パラメータ3個なので、15ビットの全件探索をすればよい。
シーケンシャルジェネレータに初期値0x7fffをセットして探索する。

4. XORSHIFT の最良のパラメータについては、Panneton and L'Ecuyer
による研究がある[2]。その研究では、シフト量だけでなく右シフトか左シフトか
などより多くの要素をパラメータにしている。
xorshift-4.cpp では 128ビット版についてのみだが、同様のパラメータの
中で最良のものを探索する。もちろん、最大周期についても確認する。
出力をソートすれば最良のパラメータが分かる。文献[2]と一致するので、
MTToolBox の計算の正当性が確認された。

5. 実装
MTToolBox はk(v)を求めたりパラメータを探索するためのものであり、
実際に疑似乱数生成器として使用するコードは別に書くことが望ましい。
高速性や可搬性に留意して求めたパラメータに対応するコードを書く
べきである。
xorshift.h と xorshift.c は 4. で求めたパラメータを使った
C 言語のソースコードである。

参考文献
[1] Marsaglia, G. 2003. Xorshift RNGs.
Journal of Statistical Software 8, 14, 1-6.
See http: //www.jstatsoft.org/v08/i14/xorshift.pdf.

[2] F. Panneton and P. L'Ecuyer,
On the Xorshift Random Number Generators,
ACM Transactions on Modeling and Computer Simulation,
2005, 15, 346-361

