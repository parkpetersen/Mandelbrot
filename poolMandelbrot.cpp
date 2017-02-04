#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <cmath>
#include <thread>
#include <mutex>
#include <queue>
#include <condition_variable>
#include <atomic>
#include <future>
#include <functional>
#include <stdexcept>
#include <memory>

struct Color
{
	int r;
	int g;
	int b;
};

template<typename T>
class TSQ
{
private:
        std::mutex m;
        std::queue<T> q;

public:
	void enqueue(T t)
	{
		std::lock_guard<std::mutex> l(m);
		q.push(t);
	}

	bool try_dequeue(T &f)
	{
		std::lock_guard<std::mutex> l(m);
		if (q.empty())
			return false;
		auto res = q.front();
		f = res;
		q.pop();
		return true;
	}

	bool empty()
    {
        if (q.empty())
        {
            return true;
        }
        return false;
    }


};

class ThreadPool2 {
    
    using func = std::function<void(void)>;
    
public:
    ThreadPool2(int);
    
    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args)
    -> std::future<typename std::result_of<F(Args...)>::type>;
    
    ~ThreadPool2();
    
private:
    
    std::vector< std::thread > pool;
    TSQ<func> queue;
    

    std::mutex itemMutex;
    std::condition_variable condition;
    std::atomic<bool> shouldContinue;

};

// the constructor just launches some amount of workers
inline ThreadPool2::ThreadPool2(int threads)
:   shouldContinue(false)
{
    for(int i = 0;i<threads;++i)
        pool.emplace_back(
                             [this]
                             {
                                 for(;;)
                                 {
                                     std::function<void(void)> task;
                                     
                                     {
                                         std::unique_lock<std::mutex> lock(this->itemMutex);
                                         this->condition.wait(lock,
                                                              [this]{ return this->shouldContinue || !this->queue.empty(); });
                                         if(this->shouldContinue && this->queue.empty())
                                             return;
                                         //task = std::move(this->tasks.front());
                                         this->queue.try_dequeue(task);
                                     }
                                     
                                     task();
                                 }
                             }
                             );
}

// add new work item to the pool
template<class F, class... Args>
auto ThreadPool2::enqueue(F&& f, Args&&... args)
-> std::future<typename std::result_of<F(Args...)>::type>
{
    using return_type = typename std::result_of<F(Args...)>::type;
    
    auto task = std::make_shared< std::packaged_task<return_type()> >(
                                                                      std::bind(std::forward<F>(f), std::forward<Args>(args)...)
                                                                      );
    
    std::future<return_type> res = task->get_future();
    {
        std::unique_lock<std::mutex> lock(itemMutex);
        
        // don't allow enqueueing after stopping the pool
        if(shouldContinue)
            throw std::runtime_error("enqueue on stopped ThreadPool");
        
        queue.enqueue([task](){ (*task)(); });
    }
    condition.notify_one();
    return res;
}

// the destructor joins all threads
inline ThreadPool2::~ThreadPool2()
{
    {
        std::unique_lock<std::mutex> lock(itemMutex);
        shouldContinue = true;
    }
    condition.notify_all();
    for(std::thread &worker: pool)
        worker.join();
}




//global variables that determine the dimensions and attributes of the images
const int imageWidth = 512;
const int imageHeight = 512;
int maxN = 255;
double minReal = -1.5;
double maxReal = 1.0;
double minImaginary = -1.0;
double maxImaginary = 1.0;
std::vector<Color> colorVect (imageHeight*imageWidth);
std::vector<std::vector<Color>> multiImageVector;
const int numThreads = 8;

//we want x to be the real number and in the mandelbrot set x scale
//this function will convert the x coordinate of each pixel
double toReal(int x, int imageWidth, double minReal, double maxReal)
{
	double range = maxReal - minReal;
	return x * (range / imageWidth) + minReal;
}

//we want y to be the imaginary number in the mandelbrot set y scale
//this function will convert the y coordinate of each pixel to imaginary
double toImaginary(int y, int imageHeight, double minImaginary, double maxImaginary)
{
	double range = maxImaginary - minImaginary;
	return y * (range / imageHeight) + minImaginary;
}

//mostly follows the wikipedia pseudocode to find how many iterations
//it takes to escape to infinity(if it does) (tests if it exceeds 2)
int findMandelbrot(double cReal, double cImaginary, int max_iterations)
{
	int i = 0;
	double zr = 0.0, zi = 0.0;
	while (i < max_iterations && zr * zr + zi * zi < 4.0)
	{
		double temp = zr * zr - zi * zi + cReal;
		zi = 2.0 * zr * zi + cImaginary;
		zr = temp;
		i++;
	}
	
	return i;
}

//assigns a color to a pixel
Color assignColor(int n)
{
	int r = (n % 256);
	int g = (n % 256);
	int b = (n % 256);
	Color *newColor = new Color();
	newColor->r = r;
	newColor->g = g;
	newColor->b = b;

	//std::cout<< newColor->r << " " << newColor->g << " " << newColor->b <<std::endl;

	return *newColor;

}

//writes the image to the file
void write()
{
	std::ofstream fout("poolMandelbrot.ppm");
	fout << "P3" << std::endl;
	fout << imageWidth << " " << imageHeight << std::endl;
	fout << "256" << std::endl;

	for (std::size_t i = 0; i < colorVect.size(); i++)
	{
		fout << colorVect[i].r << " " << colorVect[i].g << " " << colorVect[i].b << std::endl;
	}

	fout << std::endl;
	fout.close();


}

//generic template to time a function
template<typename F>
auto timeFunc(F f)
{
	auto start = std::chrono::steady_clock::now();
	f();
	auto end = std::chrono::steady_clock::now();
	return end - start;
}

//calculates and outputs the average and standard deviation
void averageStdDev(std::vector<std::chrono::duration<double>> list)
{
	int n = list.size();
	double sum = 0.0;
	for (int i = 0; i < n; i++)
	{
		sum += list[i].count();
	}
	double average = sum / n;
	std::cout << "Average: " << (sum / n) * 1000 << " milliseconds." << std::endl;
	double stdDev = 0.0;
	for (int j = 0; j < n; j++)
	{
		stdDev += pow(list[j].count() - average, 2);
	}
	std::cout << "Standard Deviation: " << sqrt(stdDev / n) * 1000 << " milliseconds." << std::endl;
}


void generateImageSection(int startingY, int endingY, int startingX, int endingX)
{
    for (int y = startingY; y < endingY; y++)
	{
		for (int x = startingX; x < endingX; x++)
		{
			double cReal = toReal(x, imageWidth, minReal, maxReal);
			double cImaginary = toImaginary(y, imageHeight, minImaginary, maxImaginary);

			int n = findMandelbrot(cReal, cImaginary, maxN);
			//std::cout<<"n: "<<n<<std::endl;
			Color pixelColor = assignColor(n);
			colorVect.at(y*512+x) = pixelColor;
		}
	}
}


//generates the image, storing it in a vector of colors. It then puts the vector into a global vector that holds all of the images being created.
void generateImageByEvenSections()
{	
	ThreadPool2 tPool(numThreads);
	std::thread threadVector[numThreads];
	for (int i = 0; i < 8; i++)
	{
		std::function<void(void)> task = std::bind(generateImageSection,(imageHeight/8)*i, (imageHeight/8)*( i + 1), 0, imageWidth);
		tPool.enqueue(task);
	} 
}

void generateImageByPixel()
{
        ThreadPool2 tPool(numThreads);
        std::thread threadVector[numThreads];
        for (int y = 0; y <imageHeight; y++)
        {
		for (int x =0; x < imageWidth ; x++)
		{
                	std::function<void(void)> task = std::bind(generateImageSection,y, y+1, x, x+1);
                	tPool.enqueue(task);
		}
        }
}

void generateImageByHalfRows()
{
        ThreadPool2 tPool(numThreads);
        std::thread threadVector[numThreads];
        for (int y = 0; y <imageHeight; y++)
        {
            	std::function<void(void)> task = std::bind(generateImageSection,y, y+1, 0 , imageWidth/2);
		tPool.enqueue(task);
		task = std::bind(generateImageSection, y, y+1, imageWidth/2, imageWidth);
		tPool.enqueue(task);
		
        }


}

void generateImageByFullRows()
{
        ThreadPool2 tPool(numThreads);
        std::thread threadVector[numThreads];
        for (int y = 0; y <imageHeight; y++)
        {
                std::function<void(void)> task = std::bind(generateImageSection,y, y+1, 0 , imageWidth);
                tPool.enqueue(task);
        }


}



int main() {

	int numImages = 5;
	std::vector<std::chrono::duration<double>> timeVector;
	std::cout << "Creating " << numImages << " Mandelbrot images." << std::endl;
	for (int i = 0; i < numImages; i++)
	{
		auto timeElapsed = timeFunc(generateImageByEvenSections);
		std::cout << "Time for image " << i+1 << ": " << std::chrono::duration_cast<std::chrono::milliseconds>(timeElapsed).count() << " milliseconds" << std::endl;
		timeVector.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(timeElapsed));
	}

	averageStdDev(timeVector);

	std::cout << "Writing images...";
	//for (int i = 0; i < numImages; i++)
	//{
		write();
	//}
	std::cout << "Done" << std::endl;
	

	std::cout << "Press enter to exit." << std::endl;
	std::cin.get();
	return 0;
}

