#pragma once

//Implementation of Shared Pointer
class ReferenceCounter;

template <typename T>
class SmartPointer
{
public:
	//default constructor
	SmartPointer() :
		pData(nullptr),
		referenceCnter(nullptr)
	{
		std::cout << "default constructor";
		std::cin.get();
		referenceCnter = new ReferenceCounter;
		referenceCnter->AddRef();
	}

	//constructor
	SmartPointer(T* pValue) :
		pData(pValue),
		referenceCnter(nullptr)
	{
		std::cout << "normal constructor";
		std::cin.get();
		referenceCnter = new ReferenceCounter;
		referenceCnter->AddRef();
	}

	//copy constructor
	SmartPointer(const SmartPointer<T>& other) :
		pData(other.pData),
		referenceCnter(other.referenceCnter)
	{
		std::cout << "copy constructor";
		std::cin.get();
		referenceCnter->AddRef();
	}

	//move constructor
	SmartPointer(SmartPointer<T>&& other)
	{
		std::cout << "move constructor";
		std::cin.get();
		pData = other.pData;
		referenceCnter = other.referenceCnter;
		referenceCnter->AddRef();

		//delete old
		other.pData = nullptr;
		other.referenceCnter->Release();
	}

	//copy assignment
	SmartPointer<T>& operator= (const SmartPointer<T>& other)
	{
		std::cout << "copy assignment";
		std::cin.get();

		if (this != &other) //avoid self assignment...no need
		{
			//if this is the last pointer
			//remove current data in the existing class
			//this deleted heap memory
			if (referenceCnter->Release() == 0)
			{
				delete pData;
				delete referenceCnter;
			}

			//copy the data from other object into existing object
			pData = other.pData;
			referenceCnter = other.referenceCnter;
			referenceCnter->AddRef();
		}
		return *this;
	}

	//move assignment
	SmartPointer<T>& operator = (SmartPointer<T>&& other)
	{
		std::cout << "move assignment";
		std::cin.get();

		if (this != &other) //avoid self assignment
		{
			//if this is the last pointer
			//remove the current data in the existing class
			//this deleted the heap memory
			if (referenceCnter->Release() == 0)
			{
				delete pData;
				delete referenceCnter;
			}

			//copy the data from the other object to existing object
			pData = other.pData;
			referenceCnter = other.referenceCnter;
			referenceCnter->AddRef();

			//delete the data in the old
			//by definition the reference count will not be zero because we just added one counter
			//so we can safely dereference
			other.referenceCnter->Release();
			other.pData = nullptr;
			other.referenceCnter = nullptr;
		}
		return *this;
	}

	//destructor
	~SmartPointer()
	{
		std::cout << "destructor";
		std::cin.get();
		//dereference
		//if last pointer - delete memory
		if (referenceCnter->Release() == 0)
		{
			delete pData;
			delete referenceCnter;
		}
	}

	//* operator
	T& operator* ()
	{
		return *pData;
	}

	//-> operator
	T* operator->()
	{
		return pData;
	}


private:
	T * pData;
	ReferenceCounter* referenceCnter;
};

class ReferenceCounter
{
public:
	void AddRef();
	int Release();

private:
	int count; //the counter
};

inline void ReferenceCounter::AddRef()
{
	count++;
}

inline int ReferenceCounter::Release()
{
	return --count; //this is a write operations
}
