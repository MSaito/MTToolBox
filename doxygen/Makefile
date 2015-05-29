doc:../docs/en/index.html ../docs/index.html

../docs:
	mkdir ../docs

../docs/index.html: mainpage.jp.txt doxygen.jp.cfg ../docs
	doxygen doxygen.jp.cfg

../docs/en/index.html: mainpage.en.txt doxygen.en.cfg ../docs
	doxygen doxygen.en.cfg