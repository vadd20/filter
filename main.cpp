#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>

namespace {
    constexpr int max = SHRT_MAX + 1;
}

class Filter {
public:
    Filter(int order, int samplesNumber, const std::vector<double> &coefficients) {
        this->order = order;
        this->samplesNumber = samplesNumber;
        this->coefficients = coefficients;
    }

    void impulse(int zeros, double value) {
        samples.clear();
        samples.resize(samplesNumber, 0);
        samples[zeros] = value;

        filter();
    }

    void step(int zeroNumber, double value) {
        samples.resize(samplesNumber, value);
        for (int i = 0; i < zeroNumber; ++i) {
            samples[i] = 0;
        }
        filter();
    }

    void output() {
        int i = 0;
        for (double &element: result) {
            printf("%d. %.15f\n", i, element);
            i++;
        }
        result.clear();
    }

private:
    int order;
    int samplesNumber;
    std::vector<double> coefficients;
    std::vector<double> samples;
    std::vector<double> result;


    // квантуем коэффициенты фильтра к целым
    std::vector<short> quantizeCoefficients() {
        std::vector<short> quantizedCoefficients(order + 1); // т.к. коэффициентов фильтра 16

        for (int i = 0; i < order + 1; ++i) {
            quantizedCoefficients[i] = (short) (coefficients[i] * max);
        }

        return quantizedCoefficients;
    }

    std::vector<short> quantizeSamples() {
        std::vector<short> quantizedSamples(samplesNumber);

        for (int i = 0; i < samplesNumber; ++i) {
            quantizedSamples[i] = (short) (samples[i] * max);
        }

        return quantizedSamples;
    }

    void convertResults(std::vector<short> &filteringResult) {
        for (int i = 0; i < samplesNumber + order; ++i) {
            result.push_back((double) filteringResult[i] / max);
        }
    }

    void filter() {
        int sum;
        std::vector<short> filteringResult; // результат после фильтрации

        auto quantizedCoefficients = quantizeCoefficients();
        auto quantizedSamples = quantizeSamples();
        for (int i = 0; i < samplesNumber + order - 1; i++) {
            sum = 0;
            for (int j = 0; j < quantizedCoefficients.size(); ++j) {
                if (i - j >= 0 && i - j < quantizedCoefficients.size()) {
                    sum += quantizedSamples[i - j] * quantizedCoefficients[j];
                }
            }
            filteringResult.push_back((short) (sum / max));
        }

        convertResults(filteringResult);
    }
};

// входной вектор настраиваемой длины (импульсная или функция хевисайда - должны задавать количество нулей, а все остальные 0.8)
// квантованные значения должны быть в filterCoefficients
int main() {
    std::vector<double> filterCoefficients = {
            -0.0052490234375,
            0.0120849609375,
            0.024200439453125,
            -0.011627197265625,
            -0.06817626953125,
            -0.021331787109375,
            0.1832275390625,
            0.391693115234375,
            0.391693115234375,
            0.1832275390625,
            -0.021331787109375,
            -0.06817626953125,
            -0.011627197265625,
            0.024200439453125,
            0.0120849609375,
            -0.0052490234375
    };

    double stepCoefficient = 0.9;
    double impulseCoefficient = 1 - pow(2, -15);

    int filterOrder = 15; // порядок фильтра
    int samplesNumberCur = 40; // количество отсчетов
    Filter filter(filterOrder, samplesNumberCur, filterCoefficients);
    int zeros = 3; // количество нулей в функции Хевисайда
    filter.step(zeros, stepCoefficient);
    std::cout << "Step response:" << std::endl;
    filter.output();
    std::cout << std::endl;

    std::cout << "Impulse response:" << std::endl;
    filter.impulse(zeros, impulseCoefficient);
    filter.output();

    return 0;
}
