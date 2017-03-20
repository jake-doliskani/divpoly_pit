
#include "Util.h"
#include <chrono>
#include <NTL/FacVec.h>

using namespace std;

Util::Util() {
}

Util::~Util() {
}

long Util::getTimeMillis() {
    auto now = chrono::system_clock::now();
    auto now_ms = chrono::time_point_cast<chrono::milliseconds>(now);
    auto value = now_ms.time_since_epoch();
    long duration = value.count();
    return duration;
}

long Util::findOrder(long a, long modulus) {
    
}

