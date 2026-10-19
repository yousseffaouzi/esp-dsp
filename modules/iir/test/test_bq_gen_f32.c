// Copyright 2026 Espressif Systems (Shanghai) PTE LTD
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <math.h>
#include "unity.h"

#include "dsps_biquad.h"
#include "dsps_biquad_gen.h"

static void allpass_frequency_response(const float *coeffs, float frequency, float *magnitude, float *phase_deg)
{
    float w = 2.0f * (float)M_PI * frequency;
    float c1 = cosf(w);
    float s1 = -sinf(w);
    float c2 = cosf(2.0f * w);
    float s2 = -sinf(2.0f * w);

    float nr = coeffs[0] + coeffs[1] * c1 + coeffs[2] * c2;
    float ni = coeffs[1] * s1 + coeffs[2] * s2;
    float dr = 1.0f + coeffs[3] * c1 + coeffs[4] * c2;
    float di = coeffs[3] * s1 + coeffs[4] * s2;
    float den = dr * dr + di * di;
    float hr = (nr * dr + ni * di) / den;
    float hi = (ni * dr - nr * di) / den;

    *magnitude = sqrtf(hr * hr + hi * hi);
    *phase_deg = atan2f(hi, hr) * 180.0f / (float)M_PI;
}

TEST_CASE("dsps_biquad_gen_allpass180_f32 frequency response", "[dsps]")
{
    const float break_frequencies[] = {0.05f, 0.1f, 0.2f, 0.4f};
    const float probe_frequencies[] = {0.01f, 0.05f, 0.1f, 0.25f, 0.49f};
    const int break_count = sizeof(break_frequencies) / sizeof(break_frequencies[0]);
    const int probe_count = sizeof(probe_frequencies) / sizeof(probe_frequencies[0]);

    for (int i = 0; i < break_count; i++) {
        float coeffs[5];
        float coeffs_other_q[5];
        float magnitude;
        float phase_deg;

        TEST_ASSERT_EQUAL(ESP_OK, dsps_biquad_gen_allpass180_f32(coeffs, break_frequencies[i], 0.5f));
        TEST_ASSERT_EQUAL(ESP_OK, dsps_biquad_gen_allpass180_f32(coeffs_other_q, break_frequencies[i], 10.0f));

        for (int j = 0; j < 5; j++) {
            TEST_ASSERT_MESSAGE(fabsf(coeffs[j] - coeffs_other_q[j]) < 1e-6f,
                                "allpass180 coefficients must not depend on qFactor");
        }

        for (int j = 0; j < probe_count; j++) {
            allpass_frequency_response(coeffs, probe_frequencies[j], &magnitude, &phase_deg);
            TEST_ASSERT_MESSAGE(fabsf(magnitude - 1.0f) < 1e-4f,
                                "allpass180 magnitude must remain unity");
        }

        allpass_frequency_response(coeffs, 0.0f, &magnitude, &phase_deg);
        TEST_ASSERT_MESSAGE(fabsf(phase_deg) < 1e-3f,
                            "allpass180 phase must be 0 degrees at DC");

        allpass_frequency_response(coeffs, break_frequencies[i], &magnitude, &phase_deg);
        TEST_ASSERT_MESSAGE(fabsf(phase_deg + 90.0f) < 1e-2f,
                            "allpass180 phase must be -90 degrees at the break frequency");

        allpass_frequency_response(coeffs, 0.499f, &magnitude, &phase_deg);
        TEST_ASSERT_MESSAGE(fabsf(phase_deg + 180.0f) < 2.0f,
                            "allpass180 phase must approach -180 degrees near Nyquist");
    }

    float coeffs180[5];
    float coeffs360[5];
    float coefficient_difference = 0.0f;
    dsps_biquad_gen_allpass180_f32(coeffs180, 0.1f, 1.0f);
    dsps_biquad_gen_allpass360_f32(coeffs360, 0.1f, 1.0f);
    for (int i = 0; i < 5; i++) {
        coefficient_difference += fabsf(coeffs180[i] - coeffs360[i]);
    }
    TEST_ASSERT_MESSAGE(coefficient_difference > 1e-3f,
                        "allpass180 and allpass360 generators must not be equivalent");
}

TEST_CASE("dsps_biquad_gen_allpass180_f32 impulse response", "[dsps]")
{
    const int len = 32;
    float input[32] = {0};
    float output[32] = {0};
    float state[2] = {0};
    float coeffs[5];

    input[0] = 1.0f;
    TEST_ASSERT_EQUAL(ESP_OK, dsps_biquad_gen_allpass180_f32(coeffs, 0.1f, 1.0f));
    TEST_ASSERT_EQUAL(ESP_OK, dsps_biquad_f32(input, output, len, coeffs, state));

    float k = -coeffs[0];
    float previous_input = 0.0f;
    float previous_output = 0.0f;

    for (int i = 0; i < len; i++) {
        float reference = k * previous_output - k * input[i] + previous_input;
        TEST_ASSERT_MESSAGE(fabsf(output[i] - reference) < 1e-5f,
                            "allpass180 processor output must match the first-order reference recursion");
        previous_input = input[i];
        previous_output = reference;
    }
}
