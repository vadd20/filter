#include <iostream>
#include <vector>
#include <cmath>

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

    void impulse(double value) {
        samples.clear();
        samples.resize(samplesNumber, 0);
        samples[0] = value;

        filter();
    }

    void step(double value) {
        samples.resize(samplesNumber, value);
        filter();
    }

    void output() {
        int i = 0;
        for (auto &element: result) {
            std::cout << i << ". " << element << "\n";
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
        for (int i = 0; i < samplesNumber; ++i) {
            result.push_back((double) filteringResult[i] / max);
        }
    }

    void filter() {
        int sum;
        std::vector<short> filteringResult; // результат после фильтрации

        auto quantizedCoefficients = quantizeCoefficients();
        auto quantizedSamples = quantizeSamples();
        for (int i = 0; i < samplesNumber; ++i) {
            sum = 0;
            for (int j = 0; j < order + 1; ++j) {
                if (i - j >= 0) {
                    sum += quantizedSamples[i - j] * quantizedCoefficients[j];
                }
            }
            filteringResult.push_back((short) (sum / max));
        }

        convertResults(filteringResult);
    }
};


int main() {
    std::vector<double> filterCoefficients = {
            -0.005240352829934454698124213223309197929,
            0.012071466320048965942257623851219250355,
            0.02421434190724541107853085009082860779,
            -0.011617204054193294715524586990795796737,
            -0.068181111645009701005548663488298188895,
            -0.021338477592380677289041202016051101964,
            0.183237394106925050030199031425581779331,
            0.391680520031719880957865598247735761106,
            0.391680520031719880957865598247735761106,
            0.183237394106925050030199031425581779331,
            -0.021338477592380677289041202016051101964,
            -0.068181111645009701005548663488298188895,
            -0.011617204054193294715524586990795796737,
            0.02421434190724541107853085009082860779,
            0.012071466320048965942257623851219250355,
            -0.005240352829934454698124213223309197929
    };

    double stepCoefficient = 0.8;
    double impulseCoefficient = 1 - pow(2, -15);

    int filterOrder = 15; // порядок фильтра
    int filterSamplesNumber = 16; // количество отсчетов
    Filter filter(filterOrder, filterSamplesNumber, filterCoefficients);
    filter.step(stepCoefficient);
    std::cout << "Step response:" << std::endl;
    filter.output();
    std::cout << std::endl;

    std::cout << "Impulse response:" << std::endl;
    filter.impulse(impulseCoefficient);
    filter.output();

    return 0;
}
