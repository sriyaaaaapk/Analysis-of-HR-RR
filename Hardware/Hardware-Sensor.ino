//sends the raw values to MATLAB
#include <Wire.h>
#include "MAX30105.h"

MAX30105 particleSensor;

const long BAUD_RATE = 115200;
const int DISPLAY_INTERVAL_MS = 50;  // 50 plotted samples/second

void setup() {
  Serial.begin(BAUD_RATE);
  delay(1000);

  if (!particleSensor.begin(Wire, I2C_SPEED_FAST)) {
    // Text error is for the Serial Monitor
    Serial.println("MAX30102_NOT_FOUND");
    while (1) {
      delay(1000);
    }
  }

  byte ledBrightness = 60;
  byte sampleAverage = 4;
  byte ledMode = 2;       // Red + IR
  int sampleRate = 400;   // Sensor sampling rate
  int pulseWidth = 411;
  int adcRange = 4096;

  particleSensor.setup(
    ledBrightness,
    sampleAverage,
    ledMode,
    sampleRate,
    pulseWidth,
    adcRange
  );

  // Keep setup messages out of the numeric plot
  delay(500);
}

void loop() {
  static unsigned long lastDisplay = 0;

  if (millis() - lastDisplay >= DISPLAY_INTERVAL_MS) {
    lastDisplay = millis();

    long irValue = particleSensor.getIR();

    // Serial Plotter expects numeric values
    Serial.println(irValue);
  }
}
