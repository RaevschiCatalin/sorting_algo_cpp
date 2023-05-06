#include <iostream>
#include <vector>
#include <cstdlib> // for rand()
#include <functional>
#include <chrono>

using namespace std;
//time
namespace mynamespace
{
    // Function to measure execution time of other functions
    template <typename Func>
    double measureExecutionTime(Func func)
    {
        // Get the current time before executing the function
        auto start = chrono::high_resolution_clock::now();

        // Execute the function
        func();

        // Get the current time after executing the function
        auto end = chrono::high_resolution_clock::now();

        // Calculate the elapsed time
        chrono::duration<double> elapsedSeconds = end - start;

        // Return the elapsed time in seconds
        return elapsedSeconds.count();
    }
}
//qick
void quickSort(vector<int>& arr, int left, int right) {
    if (left >= right) return;

    int pivot = arr[(left + right) / 2];
    int i = left, j = right;
    while (i <= j) {
        while (arr[i] < pivot) i++;
        while (arr[j] > pivot) j--;
        if (i <= j) {
            swap(arr[i], arr[j]);
            i++;
            j--;
        }
    }
    if (left < j) quickSort(arr, left, j);
    if (i < right) quickSort(arr, i, right);
}
//merge
void merge(vector<int>& arr, int left, int middle, int right) {
    vector<int> temp(right - left + 1);
    int i = left, j = middle + 1, k = 0;
    while (i <= middle && j <= right) {
        if (arr[i] <= arr[j]) {
            temp[k++] = arr[i++];
        }
        else {
            temp[k++] = arr[j++];
        }
    }
    while (i <= middle) {
        temp[k++] = arr[i++];
    }
    while (j <= right) {
        temp[k++] = arr[j++];
    }
    for (int p = 0; p < k; p++) {
        arr[left + p] = temp[p];
    }
}

void mergeSort(vector<int>& arr, int left, int right) {
    if (left < right) {
        int middle = (left + right) / 2;
        mergeSort(arr, left, middle);
        mergeSort(arr, middle + 1, right);
        merge(arr, left, middle, right);
    }
}

// Heap Sort Implementation
void heapify(vector<int>& arr, int n, int i) {
    int largest = i; // Initialize largest as root
    int l = 2 * i + 1; // left = 2*i + 1
    int r = 2 * i + 2; // right = 2*i + 2

    // If left child is larger than root
    if (l < n && arr[l] > arr[largest])
        largest = l;

    // If right child is larger than largest so far
    if (r < n && arr[r] > arr[largest])
        largest = r;

    // If largest is not root
    if (largest != i) {
        swap(arr[i], arr[largest]);

        // Recursively heapify the affected sub-tree
        heapify(arr, n, largest);
    }
}
void heapSort(vector<int>& arr) {
    int n = arr.size();

    // Build heap (rearrange array)
    for (int i = n / 2 - 1; i >= 0; i--)
        heapify(arr, n, i);

    // One by one extract an element from heap
    for (int i = n - 1; i >= 0; i--) {
        // Move current root to end
        swap(arr[0], arr[i]);

        // call max heapify on the reduced heap
        heapify(arr, i, 0);
    }
}
//insertion
void insertionSort(vector<int>& arr) {
    int n = arr.size();
    for (int i = 1; i < n; i++) {
        int key = arr[i];
        int j = i - 1;
        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
            j--;
        }
        arr[j + 1] = key;
    }
}

// Bubble Sort
void bubbleSort(vector<int>& arr) {
    int n = arr.size();
    for (int i = 0; i < n - 1; i++) {
        for (int j = 0; j < n - i - 1; j++) {
            if (arr[j] > arr[j + 1]) {
                swap(arr[j], arr[j + 1]);
            }
        }
    }
}
void countingSort(vector<int>& arr) {
    int maxVal = *max_element(arr.begin(), arr.end());
    vector<int> count(maxVal + 1);

    for (int i = 0; i < arr.size(); i++) {
        count[arr[i]]++;
    }

    int j = 0;
    for (int i = 0; i <= maxVal; i++) {
        while (count[i]-- > 0) {
            arr[j++] = i;
        }
    }
}
// Selection Sort
void selectionSort(vector<int>& arr) {
    int n = arr.size();
    for (int i = 0; i < n - 1; i++) {
        int minIndex = i;
        for (int j = i + 1; j < n; j++) {
            if (arr[j] < arr[minIndex]) {
                minIndex = j;
            }
        }
        if (minIndex != i) {
            swap(arr[i], arr[minIndex]);
        }
    }
}

// Shell Sort
void shellSort(vector<int>& arr) {
    int n = arr.size();
    for (int gap = n / 2; gap > 0; gap /= 2) {
        for (int i = gap; i < n; i++) {
            int temp = arr[i];
            int j;
            for (j = i; j >= gap && arr[j - gap] > temp; j -= gap) {
                arr[j] = arr[j - gap];
            }
            arr[j] = temp;
        }
    }
}

// Radix Sort
void radixSort(vector<int>& arr) {
    int n = arr.size();
    int maxElement = *max_element(arr.begin(), arr.end());

    for (int exp = 1; maxElement / exp > 0; exp *= 10) {
        vector<int> output(n);
        vector<int> count(10, 0);

        for (int i = 0; i < n; i++) {
            count[(arr[i] / exp) % 10]++;
        }

        for (int i = 1; i < 10; i++) {
            count[i] += count[i - 1];
        }

        for (int i = n - 1; i >= 0; i--) {
            output[count[(arr[i] / exp) % 10] - 1] = arr[i];
            count[(arr[i] / exp) % 10]--;
        }

        for (int i = 0; i < n; i++) {
            arr[i] = output[i];
        }
    }
}

// function to generate random data
vector<int> generateData(int size) {
    vector<int> data(size);
    srand(time(NULL)); // seed the random number generator with the current time
    for (int i = 0; i < size; i++) {
        data[i] = rand() % 10000000; // generate random integers between 0 and 999
    }
    return data;
}

int main() {
    int n = 100000; // number of elements in the vector
    vector<int> data = generateData(n); // generate random data
    // sort the data using Quick Sort and measure the time taken
    auto quick = mynamespace::measureExecutionTime([&data]() {
        quickSort(data, 0, data.size() - 1);
    });
    auto heappy = mynamespace::measureExecutionTime([&data]() {
        heapSort(data);
    });
    auto merge = mynamespace::measureExecutionTime([&data]() {
        mergeSort(data, 0, data.size() - 1);
    });
    auto insertionSortDuration = mynamespace::measureExecutionTime([&data]() {
        insertionSort(data);
    });
    // sort the data using Bubble Sort and measure the time taken
    auto bubbleSortDuration = mynamespace::measureExecutionTime([&data]() {
        bubbleSort(data);
    });
    auto counting = mynamespace::measureExecutionTime([&data]() {
        countingSort(data);
    });
    auto radix = mynamespace::measureExecutionTime([&data]() {
        radixSort(data);
    });
    auto selection = mynamespace::measureExecutionTime([&data]() {
        selectionSort(data);
    });
    auto shell = mynamespace::measureExecutionTime([&data]() {
        shellSort(data);
    });


    cout << "\nTime taken by Quick Sort: " << quick << " seconds\n" << endl;
    cout << "Time taken by Heap Sort: " << heappy << " seconds\n" << endl;
    cout << "Time taken by Merge Sort: " << merge << " seconds\n" << endl;
    cout << "Time taken by Insertion Sort: " << insertionSortDuration << " seconds\n" << endl;
    cout << "Time taken by Bubble Sort: " << bubbleSortDuration << " seconds\n" << endl;
    cout << "Time taken by Counting Sort: " << counting << " seconds\n" << endl;
    cout << "Time taken by Radix Sort: " << radix << " seconds\n" << endl;
    cout << "Time taken by Selection Sort: " << selection << " seconds\n" << endl;
    cout << "Time taken by Shell Sort: " << shell << " seconds\n" << endl;

    return 0;
}
