# micro-ugens

Microsound UGens for SuperCollider.
(c) 2026 by Luc Doebereiner
luc.doebereiner@gmail.com

## UGens

**SubGrain** — Subsample-accurate granular synthesis with density/regularity control, morphable envelopes, and multichannel panning. Operates on live input (internal delay line) or buffers.

**Pulsar** — Pulsar synthesis with overlapping pulsarets, polynomial envelope shaping, formant control, and Kuramoto-style phase coupling for synchronizing multiple oscillators.

## Requirements

- SuperCollider source code (for plugin headers)
- CMake ≥ 3.12
- C++17 compiler

## Building

```bash
git clone https://github.com/yourname/micro-ugens.git
cd micro-ugens
mkdir build && cd build
cmake -DSC_PATH=/path/to/supercollider ..
make
make install
```

Set `SC_PATH` to your SuperCollider source directory, or set the `SC_PATH` environment variable.

### Options

- `-DSUPERNOVA=ON` — Build for supernova server

### Install locations

- macOS: `~/Library/Application Support/SuperCollider/Extensions/MicroUGens`
- Linux: `~/.local/share/SuperCollider/Extensions/MicroUGens`
- Windows: `%LOCALAPPDATA%/SuperCollider/Extensions/MicroUGens`

After installation, recompile the class library in SuperCollider.

## License

GPL-3.0
