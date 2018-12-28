#include "test_shared_pointer.h"
#include <iostream>
#include <vector>

void test_shared_pointer() 
{
	//create new pointer
	auto x = std::make_unique<std::string>("whale"); // same as new string ("whale")
	x->append("_n");

	std::string dd = "dd";
	std::cout << dd.c_str();

	std::cout << "new word: " << x->c_str();
	std::cin.get();

	//move to a new pointer
	std::unique_ptr<std::string> x2 = std::move(x);

	//old ptr is null
	bool bln_p = x == nullptr;
	std::cout << "old pter is null: " << bln_p;
	std::cin.get();

	//unique ptr can be stored in a container
	std::vector<std::unique_ptr<std::string>> v;
	v.push_back(std::make_unique<std::string>("one"));
	v.push_back(std::make_unique<std::string>("two"));
	v.push_back(std::make_unique<std::string>("three"));

	struct B {
		virtual void bar() {
			std::cout << "B::bar";
			std::cin.get();
		}
		virtual ~B() = default;
	};

	struct D :B {
		D(const std::string& name) {
			_name = name;
			std::cout << "new D: " << _name.c_str() << "\n";
		}
		~D() { std::cout << "delete D: " << _name.c_str() << "\n"; }
		void bar() override { std::cout << "D::bar " << _name.c_str() << "\n"; }
		std::string _name;
	};

	struct C :B {
		C(std::unique_ptr<std::string> name) :
			_name(std::move(name)) {
			std::cout << "new C: " << _name->c_str() << "\n";
		}
		~C() { std::cout << "delete C:" << _name->c_str() << "\n"; }
		void bar() override { std::cout << "C::bar " << _name->c_str() << "\n"; }
		std::unique_ptr<std::string> _name;
	};

	struct G :B {
		G(std::shared_ptr<std::string> name) :
			_name(name) {
			std::cout << "new G: " << _name->c_str() << "\n";
		}
		~G() { std::cout << "delete G:" << _name->c_str() << "\n"; }
		void bar() override { std::cout << "G::bar " << _name->c_str() << "\n"; }
		std::shared_ptr<std::string> _name;
	};

	{
		auto name_x = std::make_unique<std::string>("one");
		auto name_y = std::make_unique<std::string>("two");

		//name_x = name_y; impossible

		std::vector<std::unique_ptr<B>> myV;
		myV.push_back(std::make_unique<D>("one"));
		myV.push_back(std::make_unique<D>("two"));
		//move the first one
		//myV.push_back(std::move(myV[0]));
		//add C struct
		auto myC = std::make_unique<C>(std::make_unique < std::string>("three"));
		myV.emplace_back(std::move(myC));
		//myV.emplace_back(myC); //wont work
		//myV.push_back(myC);

		for (int i = 0; i < myV.size(); ++i)
		{
			myV[i]->bar();
		}
		std::cin.get();

		//now shared pointer
		auto name_g = std::make_shared < std::string>("gFour");
		auto gObj = std::make_shared<G>(name_g);

		std::vector<std::shared_ptr<B>> myG;
		myG.push_back(gObj);

	}
}