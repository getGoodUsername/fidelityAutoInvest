#include <vector>
#include <algorithm>
#include <limits>

#include <iostream>


std::vector<double> rebalance(const std::vector<double>& targetWeights, const std::vector<double>& assetValues, double totalPortfolioValue);
std::vector<double> changeInTotalPortfolioValueNeededForCurrentAssetValuesToBecomeTheTargetValues(const std::vector<double>& targetWeights, const std::vector<double>& assetValues, double totalPortfolioValue);
bool canDoFullSingleOpRebalance(const std::vector<double>& targetWeights, const std::vector<double>& assetValues, double totalPortfolioValue, double targetBuySell);
std::size_t getMostOverWeightIndex(const std::vector<double>& targetWeights, const std::vector<double>& currAssetValues, double currPortfolioValue);
std::size_t getMostUnderWeightIndex(const std::vector<double>& targetWeights, const std::vector<double>& currAssetValues, double currPortfolioValue);
std::vector<double> singleOperationTransaction(const std::vector<double>& targetWeights, const std::vector<double>& assetValues, double targetBuySell, std::size_t maxIter = 100000);


template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
    os << "[ ";
    for (std::size_t i = 0; i < vec.size() - 1; i += 1)
    {
        os << vec[i] << ", ";
    }

    os << vec[vec.size() - 1] << "]";
    return os;
}


int main(void)
{
    std::vector<double> targetWeights = {
        36.737 / 100,
        23.369 / 100,
        8.369 / 100,
        4.021 / 100,
        6.631 / 100,
        3.631 / 100,
        17.242 / 100
    };

    std::vector<double> assetValues = {
        23672.72,
        15918.04,
        6708.47,
        3659.90,
        5608.27,
        3543.19,
        11317.59
    };

    const double targetBuySell = 20000;

    std::cout << singleOperationTransaction(targetWeights, assetValues, targetBuySell, targetBuySell) << std::endl;
    return 0;
}





std::vector<double> changeInTotalPortfolioValueNeededForCurrentAssetValuesToBecomeTheTargetValues(const std::vector<double>& targetWeights, const std::vector<double>& assetValues, double totalPortfolioValue)
{
    std::vector<double> result(targetWeights.size());
    for (std::size_t i = 0; i < targetWeights.size(); i += 1)
    {
        const double weight = targetWeights[i];
        const double value = assetValues[i];
        result[i] = ([&]() -> double
        {
            if (weight == 0) return (value > 0) ? std::numeric_limits<double>::max() : 0;
            return (value / weight) - totalPortfolioValue;
        })();
    }

    return result;
}





bool canDoFullSingleOpRebalance(const std::vector<double>& targetWeights, const std::vector<double>& assetValues, double totalPortfolioValue, double targetBuySell)
{
    const bool isSell = targetBuySell < 0;
    const std::vector<double> change = changeInTotalPortfolioValueNeededForCurrentAssetValuesToBecomeTheTargetValues(targetWeights, assetValues, totalPortfolioValue);

    if (isSell)
    {
        const double minPortfolioSellForFullOnlySellRebalance = *std::min_element(change.begin(), change.end());
        return targetBuySell <= minPortfolioSellForFullOnlySellRebalance;
    }

    const double minPortfolioBuyForFullOnlyBuyRebalance = *std::max_element(change.begin(), change.end());
    return targetBuySell >= minPortfolioBuyForFullOnlyBuyRebalance;
}





std::vector<double> rebalance(const std::vector<double>& targetWeights, const std::vector<double>& assetValues, double totalPortfolioValue)
{
    std::vector<double> result(targetWeights.size());

    for (std::size_t i = 0; i < targetWeights.size(); i += 1)
    {
        const double idealValue = totalPortfolioValue * targetWeights[i];
        result[i] = idealValue - assetValues[i];
    }

    return result;
}





std::size_t getMostOverWeightIndex(const std::vector<double>& targetWeights, const std::vector<double>& currAssetValues, double currPortfolioValue)
{
    const std::vector<double> currRebalance = rebalance(targetWeights, currAssetValues, currPortfolioValue);

    std::size_t overWeightIndex = 0;
    for (std::size_t i = 1; i < currRebalance.size(); i += 1)
    {
        overWeightIndex = (currRebalance[i] < currRebalance[overWeightIndex]) ? i : overWeightIndex;
    }
    return overWeightIndex;
}





std::size_t getMostUnderWeightIndex(const std::vector<double>& targetWeights, const std::vector<double>& currAssetValues, double currPortfolioValue)
{
    const std::vector<double> currRebalance = rebalance(targetWeights, currAssetValues, currPortfolioValue);

    std::size_t underWeightIndex = 0;
    for (std::size_t i = 1; i < currRebalance.size(); i += 1)
    {
        underWeightIndex = (currRebalance[i] > currRebalance[underWeightIndex]) ? i : underWeightIndex;
    }
    return underWeightIndex;
}





std::vector<double> singleOperationTransaction(const std::vector<double>& targetWeights, const std::vector<double>& assetValues, double targetBuySell, std::size_t maxIter)
{
    const double totalPortfolioValue = ([&]() -> double
    {
        double sum = 0;
        for (std::size_t i = 0; i < assetValues.size(); i += 1)
        {
            sum += assetValues[i];
        }

        return sum;
    })();

    if (canDoFullSingleOpRebalance(targetWeights, assetValues, totalPortfolioValue, targetBuySell))
        return rebalance(targetWeights, assetValues, totalPortfolioValue);

    const std::size_t numberOfZeroWeightAssets = ([&]()
    {
        std::size_t count = 0;
        for (const double weight : targetWeights)
        {
            count += weight == 0;
        }

        return count;
    })();


    const bool isSell = targetBuySell < 0;
    std::vector<double> endAssetValues = assetValues;
    const double changePacket = targetBuySell / maxIter;
    const auto assetSelector = (isSell) ? getMostOverWeightIndex : getMostUnderWeightIndex;
    const std::size_t initTargetAssetIndex = assetSelector(targetWeights, assetValues, totalPortfolioValue);
    const double initChange = ([&]() -> double
    {
        const double weight0 = targetWeights[initTargetAssetIndex];
        const double value0 = assetValues[initTargetAssetIndex];
        std::vector<double> assetValueChangesForEqRebalanceVal(targetWeights.size());
        const auto getValue = [&](std::size_t i) -> double
        {
            const double weight1 = targetWeights[i];
            const double value1 = assetValues[i];
            return (totalPortfolioValue * (weight1 - weight0) + value0 - value1) / (weight0 - weight1 - 1);
        };

        for (std::size_t i = 0; i < initTargetAssetIndex; i += 1)
        {
            assetValueChangesForEqRebalanceVal[i] = getValue(i);
        }

        // do this for convenience, actual value would be 0 but thats
        // an irrelevant result since want to find the value such that
        // it would make this most under/over weight index tied with
        // in place with another index.
        assetValueChangesForEqRebalanceVal[initTargetAssetIndex] = targetBuySell;

        for (std::size_t i = initTargetAssetIndex + 1; i < targetWeights.size(); i += 1)
        {
            assetValueChangesForEqRebalanceVal[i] = getValue(i);
        }

        return (isSell) ?
            *std::max_element(assetValueChangesForEqRebalanceVal.begin(), assetValueChangesForEqRebalanceVal.end()) :
            *std::min_element(assetValueChangesForEqRebalanceVal.begin(), assetValueChangesForEqRebalanceVal.end());
    })();
    double currPortfolioValue = totalPortfolioValue + initChange;
    const std::size_t remainingIter = (targetBuySell - initChange) / changePacket;

    endAssetValues[initTargetAssetIndex] += initChange;
    for (std::size_t i = 0; i < remainingIter; i += 1, currPortfolioValue += changePacket)
    {
        endAssetValues[assetSelector(targetWeights, endAssetValues, currPortfolioValue)] += changePacket;
    }

    const double targetFinalPortfolioValue = totalPortfolioValue + targetBuySell;
    endAssetValues[assetSelector(targetWeights, endAssetValues, currPortfolioValue)] += targetFinalPortfolioValue - currPortfolioValue;

    return ([&]() -> std::vector<double>
    {
        std::vector<double>& diff = endAssetValues;
        for (std::size_t i = 0; i < diff.size(); i += 1)
        {
            diff[i] -= assetValues[i];
        }

        return diff;
    })();
}