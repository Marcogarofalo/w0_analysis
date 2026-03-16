#pragma once

std::pair<std::vector<size_t>, std::vector<size_t>>
matching_indices(const std::vector<std::string>& a,
    const std::vector<std::string>& b);

    struct Match {
    std::string value;
    size_t idx1;
    size_t idx2;
    size_t idx3;
};

std::vector<Match> findAllIndices(const std::vector<std::string>& v1, 
                                  const std::vector<std::string>& v2, 
                                  const std::vector<std::string>& v3);

bool hasDuplicates(const std::vector<std::string>& vec);