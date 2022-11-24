CPP = g++ -std=c++17
CPPFLAGS = -Wall -Werror -Wextra
PROJECTNAME = s21_matrix_oop
SOURCE = $(wildcard s21_*.cpp)

all: clean $(PROJECTNAME).a test

clean:
	@rm -rf *.out
	@rm -rf *.gcda
	@rm -rf *.gcno
	@rm -rf *.a
	@rm -rf *.info
	@rm -rf *.o
	@rm -rf s21_matrix_oop.dSYM
	@rm -rf test
	@rm -rf report

$(PROJECTNAME).a: $(SOURCE)
	@$(CPP) $(CPPFLAGS) $(SOURCE) -c
	@ar rcs $(PROJECTNAME).a *.o
	@ranlib $(PROJECTNAME).a

test:
	@clear
	@$(CPP) $(CPPFLAGS) $(PROJECTNAME)_test.cpp $(PROJECTNAME).a -o test -lgtest
	./test

gcov_report:
	@$(CPP) $(CPPFLAGS) -fprofile-arcs -ftest-coverage $(PROJECTNAME)_test.cpp $(PROJECTNAME).cpp -o test -lgtest
	@./test
	@lcov -t "test" -o test.info -c -d .
	@genhtml -o report test.info
	@open ./report/index.html

check:
	clang-format -n *.cpp
	clang-format -n *.h