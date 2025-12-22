# Concepts for RADAR LAB Exam

---

# 1️⃣ Continuous Wave (CW) Radar

![Image](https://www.tutorialspoint.com/radar_systems/images/cw_radar.jpg?utm_source=chatgpt.com)

![Image](https://upload.wikimedia.org/wikipedia/commons/3/30/Cw_radar.png?utm_source=chatgpt.com)

![Image](https://ars.els-cdn.com/content/image/3-s2.0-B9780081004050000136-f08-02-9780081004050.jpg?utm_source=chatgpt.com)

## Concept

A **CW radar** transmits a **continuous sinusoidal signal** and listens for its echo.
It **cannot measure range**, only **velocity**, because there is no time reference.

## Principle

If the target is moving, the reflected signal experiences a **Doppler shift**.

### Doppler Frequency

$$
f_d = \frac{2v}{\lambda}
$$

where

* (v) = target velocity
* $$(\lambda = \frac{c}{f_c})$$

## Key Points for Exam

* Measures **velocity only**
* Simple hardware
* Used in **speed guns, motion sensors**

### MATLAB Example (Velocity Detection)

```matlab
fc = 10e9;      % Carrier frequency (10 GHz)
c = 3e8;
v = 30;         % Target speed (m/s)

lambda = c/fc;
fd = 2*v/lambda;

t = 0:1e-6:1e-3;
rx = cos(2*pi*fd*t);

plot(t, rx)
xlabel('Time (s)')
ylabel('Amplitude')
title('Doppler Shift in CW Radar')
```

---

# 2️⃣ Pulse Radar

## Concept

Pulse radar transmits **short pulses** and measures **echo delay** to compute **range**.

## Range Equation

$$
R = \frac{c \cdot T}{2}
$$

where (T) is round-trip delay.

## Advantages

* Measures **range**
* Can measure **velocity** using Doppler processing

### MATLAB Example (Range Measurement)

```matlab
c = 3e8;
R = 5000; % meters
delay = 2*R/c;

fs = 1e6;
t = 0:1/fs:0.001;
tx = pulstran(t,0, @rectpuls,1e-5);
rx = pulstran(t,delay, @rectpuls,1e-5);

plot(t, tx, t, rx)
legend('Tx','Rx')
title('Pulse Radar Echo Delay')
```

---

# 3️⃣ Velocity (Doppler Radar)

![Image](https://electronicsdesk.com/wp-content/uploads/2019/08/geometry-of-radar-and-target-in-deriving-doppler-frequency-shift.jpg?utm_source=chatgpt.com)

![Image](https://amt.copernicus.org/articles/15/7315/2022/amt-15-7315-2022-f01.png?utm_source=chatgpt.com)

## Concept

Velocity causes **frequency shift** in reflected signal.

## Formula

$$
v = \frac{f_d \lambda}{2}
$$

### MATLAB Example

```matlab
fd = 200;        % Doppler frequency
fc = 10e9;
lambda = 3e8/fc;

v = fd*lambda/2
```

---

# 4️⃣ Range

## Concept

Distance between radar and target.

## Depends On

* Pulse width
* Sampling frequency

### Range Resolution

$$
\Delta R = \frac{c \tau}{2}
$$

### MATLAB Example

```matlab
tau = 1e-6;
c = 3e8;

range_resolution = c*tau/2
```

---

# 5️⃣ Azimuth

## Concept

Azimuth is **horizontal angle** of target.

## Measured Using

* Antenna rotation
* Phased array beam steering

### MATLAB (Angle Estimation – Simplified)

```matlab
theta = -60:0.1:60;
beam = cosd(theta).^4;

plot(theta, beam)
xlabel('Azimuth Angle (deg)')
ylabel('Beam Pattern')
title('Azimuth Beam Pattern')
```

---

# 6️⃣ Elevation

## Concept

Vertical angle of target.

## Similar to Azimuth

But measured in **vertical plane**.

### MATLAB Example

```matlab
theta = -30:0.1:30;
pattern = cosd(theta).^6;

plot(theta, pattern)
xlabel('Elevation Angle (deg)')
ylabel('Gain')
title('Elevation Pattern')
```

---

# 7️⃣ Matched Filter

![Image](https://www.radartutorial.eu/10.processing/pic/matched_filter.print.png?utm_source=chatgpt.com)

![Image](https://wirelesspi.com/wp-content/uploads/2024/06/figure-pulse-compression-effect-on-matched-filter-output-when-objects-are-closely-spaced.gif?utm_source=chatgpt.com)

![Image](https://i.sstatic.net/XJrJp.png?utm_source=chatgpt.com)

## Concept

A **matched filter maximizes SNR** in noisy radar signals.

## Key Idea

Filter impulse response is **time-reversed conjugate** of transmitted signal.

$$
h(t) = s^*(T - t)
$$

### MATLAB Example

```matlab
t = 0:1e-6:50e-6;
tx = chirp(t,1e6,t(end),5e6);
rx = tx + 0.5*randn(size(tx));

mf = fliplr(conj(tx));
out = conv(rx, mf);

plot(abs(out))
title('Matched Filter Output')
```

---

# 8️⃣ Kalman Filter (Radar Tracking)

![Image](https://www.skyradar.com/hs-fs/hubfs/Images/Products/FreeScopes-ATC-II/figure-2-block-diagram-of-tracking-feature.png?name=figure-2-block-diagram-of-tracking-feature.png\&width=1062\&utm_source=chatgpt.com)

## Concept

Used to **track target position and velocity** in noise.

## State Vector

$$
x = [position; velocity]^T
$$

### MATLAB Example

```matlab
A = [1 1; 0 1];
H = [1 0];
Q = eye(2)*0.01;
R = 1;

x = [0;1];
P = eye(2);

z = (1:50) + randn(1,50);

for k = 1:50
    x = A*x;
    P = A*P*A' + Q;

    K = P*H'/(H*P*H' + R);
    x = x + K*(z(k) - H*x);
    P = (eye(2)-K*H)*P;

    est(k) = x(1);
end

plot(z,'r.'); hold on;
plot(est,'b')
legend('Measurement','Kalman Estimate')
```

---

# 9️⃣ Noise

## Concept

Unwanted random signal affecting detection.

## Types

* Thermal noise
* Receiver noise

### MATLAB Example

```matlab
t = 0:0.001:1;
signal = sin(2*pi*5*t);
noise = randn(size(t))*0.3;

plot(t, signal+noise)
title('Signal with Noise')
```

---

# 🔟 Clutter

![Image](https://upload.wikimedia.org/wikipedia/commons/7/70/Radar-artefacts2.PNG?utm_source=chatgpt.com)

![Image](https://radarscope.zendesk.com/hc/article_attachments/4416723055122/Ground_Clutter_3.png?utm_source=chatgpt.com)

## Concept

Unwanted echoes from **ground, sea, buildings**.

## Characteristics

* Low Doppler
* High power

### MATLAB Example

```matlab
clutter = 5*rand(1,1000);
target = zeros(1,1000);
target(500) = 20;

plot(clutter + target)
title('Radar Return with Clutter')
```

---

# 1️⃣1️⃣ MIMO Radar

![Image](https://www.nwengineeringllc.com/article/img/mimo-radar-1.png?utm_source=chatgpt.com)

![Image](https://www.mathworks.com/help/examples/phased/win64/WaveformDesignDualFunctionMIMORadComSystemExample_07.png?utm_source=chatgpt.com)

![Image](https://pub.mdpi-res.com/remotesensing/remotesensing-16-01016/article_deploy/html/images/remotesensing-16-01016-g003.png?1711332382=\&utm_source=chatgpt.com)

## Concept

Multiple transmit and receive antennas → **better resolution & diversity**.

## Benefits

* Improved angle estimation
* Higher robustness

### MATLAB Example (Conceptual)

```matlab
tx = [1 0; 0 1]; % Orthogonal waveforms
rx = tx + randn(2);

imagesc(abs(rx))
colorbar
title('MIMO Radar Signal Matrix')
```

---

## 🧠 Exam Tips (Very Important)

* **CW radar → velocity only**
* **Pulse radar → range**
* **Matched filter → SNR maximization**
* **Kalman filter → tracking**
* **Clutter ≠ noise**
* **MIMO → spatial diversity**

---
