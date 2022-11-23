CPP = g++ -std=c++17
CPPFLAGS = -Wall -Werror -Wextra
PROJECTNAME = s21_matrix_oop
SOURCE = $(wildcard s21_*.cpp)

all: clean s21_matrix_oop.a test

clean:
	@rm -rf *.out
	@rm -rf $(PROJECTNAME).a
	@rm -rf $(PROJECTNAME).o
	@rm -rf s21_matrix_oop.dSYM
	@rm -rf test

s21_matrix_oop.a: $(SOURCE)
	@$(CPP) $(CPPFLAGS) $(SOURCE) -c
	@ar rcs $(PROJECTNAME).a *.o
	@ranlib $(PROJECTNAME).a

test:
	@clear
	@$(CPP) $(CPPFLAGS) test.cpp $(PROJECTNAME).a -o test  -lgtest 
	@leaks -q -atExit -- ./test

check:
	@clang-format -n -i *.cpp
	@clang-format -n -i *.h

push:
	git add *.cpp
	git add *.h
	git add .clang-format
	git commit -m "update"
	git push origin develop