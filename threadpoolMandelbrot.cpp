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
		std::cout<<"hi3"<<std::endl;

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


};

class ThreadPool
{
using func = std::function<void(void)>;
private:
        std::vector<std::thread> pool;
        TSQ<func> queue;
        std::condition_variable hasItem;
        std::mutex itemMutex;
        std::atomic<bool> shouldContinue;

public:
	ThreadPool(int n):pool(n), queue(), hasItem(), itemMutex(), shouldContinue(true)
	{
		//for(auto&t:pool)
	//	{
	//		t = std::thread([=]{this->run();});
	//	}
	}
	
	void stop(){shouldContinue = false; hasItem.notify_all();}
	void post(func f){queue.enqueue(f); hasItem.notify_one();}
	void run()
	{
		while(shouldContinue)
		{
			func f;
			while(!queue.try_dequeue(f))
			{
				std::cout<<"doing that thing"<<std::endl;
				std::unique_lock<std::mutex> l(itemMutex);
				hasItem.wait(l);
				if(!shouldContinue) return;
			}
			std::cout<<"hi"<<std::endl;
			f();
		}
	}
	void start()
	{
		for(auto&t:pool)
                {
                        t = std::thread([=]{this->run();});
                }

	}


};
//global variables that determine the dimensions and attributes of the images
const int imageWidth = 512;
const int imageHeight = 512;
int maxN = 255;
double minReal = -1.5;
double maxReal = 1.0;
double minImaginary = -1.0;
double maxImaginary = 1.0;
std::vector<std::vector<Color>> multiImageVector;

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
void write(std::vector<Color> colorVect, const int imageWidth, const int imageHeight, int i)
{
	std::ofstream fout("poolMandelbrot" +std::to_string(i) +".ppm");
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


void generateImageSection(std::vector<Color> &colorVect, int startingY, int endingY)
{
std::cout<<"hi no 2"<<std::endl;
    for (int y = startingY; y < endingY; y++)
	{
		for (int x = 0; x < imageWidth; x++)
		{
			double cReal = toReal(x, imageWidth, minReal, maxReal);
			double cImaginary = toImaginary(y, imageHeight, minImaginary, maxImaginary);

			int n = findMandelbrot(cReal, cImaginary, maxN);
			std::cout<<"n: "<<n<<std::endl;
			Color pixelColor = assignColor(n);
			colorVect.at(y*512+x) = pixelColor;
		}
	}
}


//generates the image, storing it in a vector of colors. It then puts the vector into a global vector that holds all of the images being created.
void generateImage()
{	
	int numThreads = 8;
	ThreadPool tPool(numThreads);
	std::vector<Color> colorVect (imageHeight*imageWidth);
	std::thread threadVector[numThreads];
	for (int i = 0; i < numThreads; i++)
	{
		std::function<void(void)> task = std::bind(generateImageSection, std::ref(colorVect), (imageHeight/8)*i, (imageHeight/8)*( i + 1));
		tPool.post(task);
	} 
	tPool.start();

	for (int i =0; i < imageHeight*imageWidth; i++)
	{
	//	threadVector[i].join();
	}	
	multiImageVector.push_back(colorVect);
}

int main() {

	int numImages = 5;
	std::vector<std::chrono::duration<double>> timeVector;
	std::cout << "Creating " << numImages << " Mandelbrot images." << std::endl;
	for (int i = 0; i < numImages; i++)
	{
		auto timeElapsed = timeFunc(generateImage);
		std::cout << "Time for image " << i+1 << ": " << std::chrono::duration_cast<std::chrono::milliseconds>(timeElapsed).count() << " milliseconds" << std::endl;
		timeVector.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(timeElapsed));
	}

	averageStdDev(timeVector);

	std::cout << "Writing images...";
	for (int i = 0; i < numImages; i++)
	{
		write(multiImageVector[i], imageWidth, imageHeight, i);
	}
	std::cout << "Done" << std::endl;
	

	std::cout << "Press enter to exit." << std::endl;
	std::cin.get();
	return 0;
}

