SIMD oriented Fast Mersenne Twister sample
==========================================

sfmtdc SFMTのパラメータを生成する
usage:
./sfmtdc [-s seed] [-v] [-c count] [-f outputfile] mexp

--verbose, -v        計算時間とか表示する。
--file, -f filename  探索されたパラメータはこのファイルに出力される。指定されなければ
	   	     標準出力に出る。
--count, -c count    出力件数
--seed, -s seed      乱数の種
mexp                 メルセンヌ指数

calc_equidist SFMTのパラメータから均等分布次元を計算する
usage:
./calc_equidist [-v] mexp,pos1,sl1,sl2,sr1,sr2,msk1,msk2,msk3,msk4,parity1,parity2,parity3,parity4

--verbose, -v        各v毎のvビット精度均等分布次元を出力する。
外の引数はsfmtdc の出力を使う。

test_linearity 線形性のテスト（デバッグ用）
./test_linearity
 mexp,pos1,sl1,sl2,sr1,sr2,msk1,msk2,msk3,msk4,parity1,parity2,parity3,parity4
引数はsfmtdc の出力を使う。

ちょっといい話かも知れない話
sfmt のパラメータのシフト量は、instrincs関数で書くと、リテラルにしなければならない訳だが、
searchsfmt.hpp のsetUpParametersメソッドを修正してシフト量を固定値にしてしまえば、
シフト量が同じSFMTのパラーメータを生成できる。そうなればSFMTの疑似乱数生成コードをコンパイル
しなおさなくても、（シフト量以外の）パラメータの異なるSFMTに対応できる。
ただし、固定したシフト量の値によっては、パラメータがまったく検索されないという可能性もある。
