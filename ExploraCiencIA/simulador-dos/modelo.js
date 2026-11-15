/*
 * modelo.js — Lógica de un simulador de dinámica poblacional (sin interfaz)
 * ============================================================================
 * JavaScript puro, sin dependencias, sin DOM, sin red.
 * Funciona abriendo el archivo con <script src="modelo.js"> (file://) y en Node.js.
 *
 *   Navegador:  <script src="modelo.js"></script>  → objeto global ModeloPoblacional
 *               Pruebas: ModeloPoblacional.pruebas()   (en la consola)
 *   Node:       const M = require('./modelo.js');
 *               Pruebas: node modelo.js --test
 *
 * MODELO (tiempo discreto, un paso = una unidad de tiempo)
 * --------------------------------------------------------
 *   N(t+1) = N(t) + B − D + I − E
 *
 *   B = b_ef · N      natalidad per cápita (hijos por individuo por paso, b ≥ 0)
 *   D = d · N         mortalidad per cápita (probabilidad de morir por paso, d ∈ [0,1])
 *   I = constante     inmigración ABSOLUTA (individuos que llegan por paso, I ≥ 0)
 *   E = e · N         emigración PER CÁPITA (probabilidad de emigrar por paso, e ∈ [0,1])
 *
 * DECISIONES (todas justificadas en la sección "Decisiones" más abajo):
 *   1. Emigración per cápita (no absoluta).
 *   2. Restricción d + e ≤ 1 (muerte y emigración son riesgos que compiten).
 *   3. Todos los flujos del paso se calculan con N(t); los inmigrantes de este paso
 *      no nacen, no mueren ni emigran hasta el paso siguiente.
 *   4. Dependencia de la densidad (opcional, apagada por defecto) solo sobre la natalidad:
 *        b_ef = b · max(0, 1 − N/K)      → siempre en [0, b]
 *   5. Extinción = estado terminal: si N = 0 no revive (tampoco por inmigración).
 *   6. Determinista: N real (no se redondea). Estocástico: N entero.
 *
 * ============================================================================
 * API
 * ============================================================================
 * Parámetros (objeto):
 *   modo        'determinista' (por defecto) | 'estocastico'
 *   N0          población inicial (OBLIGATORIO). Determinista: 0 o ≥ 1. Estocástico: entero ≥ 0.
 *   natalidad   b, tasa per cápita por paso, 0 ≤ b ≤ 100 (OBLIGATORIO)
 *   mortalidad  d, en [0,1] (OBLIGATORIO)
 *   inmigracion I, individuos por paso, ≥ 0 (por defecto 0)
 *   emigracion  e, en [0,1] (por defecto 0). Debe cumplirse d + e ≤ 1.
 *   densidad    true | false (por defecto false)
 *   K           capacidad de carga, > 0 (OBLIGATORIO solo si densidad = true)
 *   semilla     entero para el generador pseudoaleatorio (por defecto 1; solo se usa en estocástico)
 *
 * Funciones:
 *   validarParametros(p)  → { ok, errores[], advertencias[], parametros }
 *   crearEstado(p)        → estado inicial (lanza Error con mensaje claro si p es inválido)
 *   paso(estado)          → { estado: nuevoEstado, registro }
 *   simular(estado, n)    → { parametros, estadoFinal, trayectoria[] }  (n+1 registros: t=0..n)
 *   exportarCSV(resultado, opciones?) → string
 *   equilibrio(p)         → { N, estable } | null   (equilibrio determinista teórico)
 *   pruebas(opciones?)    → { total, pasaron, fallaron, resultados[] }
 *
 * Estado:   { t, N, extinta, parametros, rng }   (objeto congelado; nunca se muta)
 * Registro: { t, N, B, D, I, E, extinta, aviso }
 *   - t, N: tiempo y población DESPUÉS del paso.
 *   - B, D, I, E: valores efectivamente aplicados en ese paso (no los nominales).
 *     En el registro inicial (t0) valen null porque no hubo transición.
 *   - aviso: null o texto (saturación del límite numérico, extinción por umbral).
 *
 * ============================================================================
 * DECISIONES
 * ============================================================================
 * Emigración per cápita. Emigrar es una decisión de los individuos presentes: cuantos
 * más hay, más se van. Una emigración absoluta fija exigiría "recortarla" cuando N es
 * pequeño (E ≤ N) y produciría tasas per cápita enormes en poblaciones pequeñas, lo
 * que no tiene sentido biológico. Con e ∈ [0,1], E = e·N ≤ N por construcción.
 * La inmigración sí es absoluta porque depende de la fuente externa, no de N.
 *
 * Riesgos que compiten (d + e ≤ 1). Un individuo no puede morir y emigrar a la vez.
 * Si d + e > 1, D + E podría superar a N. Con d + e ≤ 1 se garantiza D + E ≤ N,
 * por lo que N(t+1) ≥ 0 siempre. En el modo estocástico se muestrea primero
 * D ~ Binomial(N, d) y luego E ~ Binomial(N − D, e/(1 − d)); así, la probabilidad
 * marginal de que un individuo emigre es exactamente e y nunca emigra un muerto.
 *
 * Natalidad estocástica: B ~ Poisson(b_ef · N). No se usa binomial porque b puede
 * ser > 1 (camadas). Es una aproximación: supone nacimientos independientes.
 * Inmigración estocástica: I ~ Poisson(I) (llegadas independientes con tasa media I).
 *
 * Muestreo eficiente con N muy grande: binomial por inversión exacta si n·p < 30
 * (tras reflejar p > 0,5); si no, aproximación normal redondeada y acotada a [0, n]
 * (error despreciable porque n·p·(1−p) ≥ 15). Poisson: método de Knuth si λ < 30; si no,
 * aproximación normal redondeada y acotada a ≥ 0. Costo por paso O(1).
 *
 * Extinción. N = 0 es terminal: la inmigración NO reactiva una población extinta.
 * Los pasos siguientes devuelven N = 0 y flujos 0. En el modo determinista N es real
 * y nunca llegaría exactamente a 0; se considera extinta cuando N < 1 individuo
 * (umbral = 1) y se fija N = 0 (el registro lleva un aviso en ese paso).
 * Como los inmigrantes del paso llegan en el mismo paso, mientras I > 0 una población
 * no se extingue en la práctica (determinista) o solo con probabilidad e^(−I) (estocástico).
 *
 * Límite numérico. N se acota a Number.MAX_SAFE_INTEGER (2^53 − 1 ≈ 9,0·10^15). Si el
 * paso lo superaría, N se fija en ese tope y el registro lleva un aviso. N0 mayor se rechaza.
 *
 * Generador: mulberry32, con estado de 32 bits guardado dentro del estado de la simulación
 * (no hay estado global). Semillas que difieran en múltiplos de 2^32 son equivalentes.
 *
 * Equilibrio determinista (para contrastar con las pruebas):
 *   Sin densidad: r = b − d − e. Si r < 0: N* = I / (−r) (estable). Si r ≥ 0: no hay (null).
 *   Con densidad: N* es la raíz positiva de (b/K)·N² − (b − d − e)·N − I = 0.
 *   Con I = 0: N* = K·(1 − (d + e)/b). Estable si |1 + r − 2·b·N* / K| < 1.
 *   El mapa discreto puede oscilar o ser caótico si r = b − d − e ≥ 2.
 * ============================================================================
 */
(function (raiz) {
  'use strict';

  // ---------------------------------------------------------------- constantes
  var LIMITE_N = Number.MAX_SAFE_INTEGER;   // tope numérico de N
  var UMBRAL_EXTINCION = 1;                 // determinista: N < 1 individuo ⇒ extinta
  var UMBRAL_NORMAL = 30;                   // por debajo: muestreo exacto; por encima: normal
  var B_MAX = 100;                          // cota de cordura para la natalidad per cápita
  var TOL_SUMA = 1e-12;                     // tolerancia de coma flotante para d + e ≤ 1
  var MODOS = ['determinista', 'estocastico'];
  var CLAVES = ['modo', 'N0', 'natalidad', 'mortalidad', 'inmigracion', 'emigracion',
                'densidad', 'K', 'semilla'];

  // ------------------------------------------------------------------ aleatorio
  // Generador mulberry32. Devuelve un objeto local; el estado de 32 bits se extrae
  // con estado() y se guarda en el estado de la simulación (sin variables globales).
  function crearRng(semilla) {
    var a = semilla >>> 0;
    return {
      siguiente: function () {               // uniforme en [0, 1)
        a = (a + 0x6D2B79F5) | 0;
        var t = Math.imul(a ^ (a >>> 15), 1 | a);
        t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
      },
      estado: function () { return a >>> 0; }
    };
  }

  // Normal estándar (Box–Muller). Descarta la segunda variable para no guardar estado extra.
  function normal(rng) {
    var u1 = 1 - rng.siguiente();            // en (0, 1]: evita log(0)
    var u2 = rng.siguiente();
    return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
  }

  // Poisson(lambda) → entero ≥ 0.
  function poisson(rng, lambda) {
    if (!(lambda > 0)) return 0;
    if (lambda < UMBRAL_NORMAL) {            // Knuth: exacto
      var L = Math.exp(-lambda), k = 0, prod = 1;
      do { k++; prod *= rng.siguiente(); } while (prod > L);
      return k - 1;
    }
    var x = Math.round(lambda + Math.sqrt(lambda) * normal(rng));  // aproximación normal
    return x < 0 ? 0 : x;
  }

  // Binomial(n, p) → entero en [0, n]. n entero ≥ 0, p en [0,1].
  function binomial(rng, n, p) {
    if (n <= 0 || p <= 0) return 0;
    if (p >= 1) return n;
    if (p > 0.5) return n - binomial(rng, n, 1 - p);    // reflejo: trabaja con p ≤ 0,5
    var q = 1 - p, media = n * p;
    if (media < UMBRAL_NORMAL) {                         // inversión exacta
      var u = rng.siguiente();
      var f = Math.exp(n * Math.log1p(-p));              // P(X = 0) = q^n
      var r = p / q, k = 0;
      while (u > f && k < n && k < 1000) {
        u -= f; k++;
        f *= ((n - k + 1) / k) * r;                      // recurrencia de la pmf
      }
      return k;
    }
    var x = Math.round(media + Math.sqrt(n * p * q) * normal(rng)); // aproximación normal
    return x < 0 ? 0 : (x > n ? n : x);
  }

  // ----------------------------------------------------------------- validación
  function esNum(x) { return typeof x === 'number' && isFinite(x); }

  function validarParametros(entrada) {
    var errores = [], advertencias = [];
    if (entrada === null || typeof entrada !== 'object' || Array.isArray(entrada)) {
      return { ok: false, errores: ['Los parámetros deben ser un objeto.'], advertencias: advertencias, parametros: null };
    }
    Object.keys(entrada).forEach(function (k) {
      if (CLAVES.indexOf(k) < 0) advertencias.push('Parámetro desconocido ignorado: "' + k + '".');
    });
    function tomar(k, defecto) { return entrada[k] === undefined ? defecto : entrada[k]; }

    var modo = tomar('modo', 'determinista');
    if (MODOS.indexOf(modo) < 0) errores.push('modo debe ser "determinista" o "estocastico" (recibido: ' + String(modo) + ').');

    var N0 = tomar('N0', undefined);
    if (N0 === undefined) errores.push('N0 (población inicial) es obligatorio.');
    else if (!esNum(N0) || N0 < 0) errores.push('N0 debe ser un número finito ≥ 0 (recibido: ' + String(N0) + ').');
    else if (N0 > LIMITE_N) errores.push('N0 supera el máximo representable con exactitud (' + LIMITE_N + ').');
    else if (modo === 'estocastico' && Math.floor(N0) !== N0) {
      advertencias.push('N0 = ' + N0 + ' no es entero; en modo estocástico se redondea a ' + Math.round(N0) + '.');
      N0 = Math.round(N0);
    } else if (modo === 'determinista' && N0 > 0 && N0 < UMBRAL_EXTINCION) {
      errores.push('En modo determinista N0 debe ser 0 o ≥ 1 individuo (recibido: ' + N0 + ').');
    }
    if (N0 === 0) advertencias.push('N0 = 0: la población está extinta desde el inicio y no se reactiva (ni por inmigración).');

    var b = tomar('natalidad', undefined);
    if (b === undefined) errores.push('natalidad es obligatoria.');
    else if (!esNum(b) || b < 0 || b > B_MAX) errores.push('natalidad debe ser un número en [0, ' + B_MAX + '] (recibido: ' + String(b) + ').');

    var d = tomar('mortalidad', undefined);
    if (d === undefined) errores.push('mortalidad es obligatoria.');
    else if (!esNum(d) || d < 0 || d > 1) errores.push('mortalidad debe estar en [0, 1] (recibido: ' + String(d) + ').');

    var I = tomar('inmigracion', 0);
    if (!esNum(I) || I < 0 || I > LIMITE_N) errores.push('inmigracion debe ser un número finito ≥ 0 (recibido: ' + String(I) + ').');

    var e = tomar('emigracion', 0);
    if (!esNum(e) || e < 0 || e > 1) errores.push('emigracion debe estar en [0, 1] (recibido: ' + String(e) + ').');

    if (esNum(d) && esNum(e) && d >= 0 && e >= 0 && d <= 1 && e <= 1 && d + e > 1 + TOL_SUMA) {
      errores.push('mortalidad + emigracion no puede superar 1 (recibido: ' + (d + e) + '): un individuo no puede morir y emigrar a la vez.');
    }

    var dens = tomar('densidad', false);
    if (typeof dens !== 'boolean') errores.push('densidad debe ser true o false (recibido: ' + String(dens) + ').');
    var K = tomar('K', null);
    if (dens === true) {
      if (K === null || K === undefined) errores.push('K (capacidad de carga) es obligatoria cuando densidad = true.');
      else if (!esNum(K) || K <= 0) errores.push('K debe ser un número finito > 0 (recibido: ' + String(K) + ').');
    } else if (K !== null && K !== undefined) {
      advertencias.push('K se ignora porque densidad = false.');
      K = null;
    }
    if (dens === true && esNum(b) && esNum(d) && esNum(e) && (b - d - e) >= 2) {
      advertencias.push('Con densidad activada y natalidad − mortalidad − emigracion ≥ 2 el modelo discreto puede oscilar o ser caótico.');
    }

    var semilla = tomar('semilla', 1);
    if (!esNum(semilla)) errores.push('semilla debe ser un número entero finito (recibido: ' + String(semilla) + ').');
    else if (Math.floor(semilla) !== semilla) {
      advertencias.push('semilla = ' + semilla + ' no es entera; se trunca a ' + Math.trunc(semilla) + '.');
      semilla = Math.trunc(semilla);
    }

    if (errores.length) return { ok: false, errores: errores, advertencias: advertencias, parametros: null };
    var normalizados = Object.freeze({
      modo: modo, N0: N0, natalidad: b, mortalidad: d, inmigracion: I, emigracion: e,
      densidad: dens, K: dens ? K : null, semilla: semilla
    });
    return { ok: true, errores: [], advertencias: advertencias, parametros: normalizados };
  }

  // ------------------------------------------------------------------- estado
  function crearEstado(parametros) {
    var v = validarParametros(parametros);
    if (!v.ok) throw new Error('Parámetros inválidos:\n - ' + v.errores.join('\n - '));
    var p = v.parametros;
    return Object.freeze({
      t: 0, N: p.N0, extinta: p.N0 === 0, parametros: p,
      rng: p.semilla >>> 0
    });
  }

  function comprobarEstado(estado) {
    if (!estado || typeof estado !== 'object' || !estado.parametros || !esNum(estado.N) ||
        !esNum(estado.t) || estado.N < 0) {
      throw new Error('Estado inválido: use crearEstado(parametros) para obtener uno.');
    }
  }

  // Un paso. Función pura: no modifica "estado"; devuelve uno nuevo y el registro.
  function paso(estado) {
    comprobarEstado(estado);
    var p = estado.parametros, t1 = estado.t + 1;

    if (estado.extinta) {                      // estado terminal
      return {
        estado: Object.freeze({ t: t1, N: 0, extinta: true, parametros: p, rng: estado.rng }),
        registro: Object.freeze({ t: t1, N: 0, B: 0, D: 0, I: 0, E: 0, extinta: true, aviso: null })
      };
    }

    var N = estado.N;
    // Natalidad efectiva (densidad): b·max(0, 1 − N/K) ∈ [0, b]
    var bEf = p.densidad ? p.natalidad * Math.max(0, 1 - N / p.K) : p.natalidad;
    var B, D, I, E, rngNuevo = estado.rng;

    if (p.modo === 'determinista') {
      B = bEf * N; D = p.mortalidad * N; E = p.emigracion * N; I = p.inmigracion;
    } else {
      var rng = crearRng(estado.rng);
      D = binomial(rng, N, p.mortalidad);
      var pe = p.mortalidad >= 1 ? 0 : Math.min(1, p.emigracion / (1 - p.mortalidad));
      E = binomial(rng, N - D, pe);
      B = poisson(rng, bEf * N);
      I = poisson(rng, p.inmigracion);
      rngNuevo = rng.estado();
    }

    var Nn = (N - D - E) + (B + I);
    var aviso = null;
    if (Nn < 0) Nn = 0;                        // solo por redondeo de coma flotante
    if (!(Nn <= LIMITE_N)) {
      Nn = LIMITE_N;
      aviso = 'N saturada en el máximo numérico (' + LIMITE_N + '); N ya no cumple N(t+1) = N + B − D + I − E.';
    }
    var extinta = false;
    if (p.modo === 'determinista' ? Nn < UMBRAL_EXTINCION : Nn === 0) {
      if (Nn > 0) aviso = 'Extinción por umbral: N < 1 individuo en el modelo determinista.';
      Nn = 0; extinta = true;
    }

    return {
      estado: Object.freeze({ t: t1, N: Nn, extinta: extinta, parametros: p, rng: rngNuevo }),
      registro: Object.freeze({ t: t1, N: Nn, B: B, D: D, I: I, E: E, extinta: extinta, aviso: aviso })
    };
  }

  // n pasos. Devuelve n+1 registros (el primero es el estado inicial, con flujos null).
  function simular(estado, n) {
    comprobarEstado(estado);
    if (!esNum(n) || n < 0 || Math.floor(n) !== n) throw new Error('n debe ser un entero ≥ 0 (recibido: ' + String(n) + ').');
    var trayectoria = [Object.freeze({
      t: estado.t, N: estado.N, B: null, D: null, I: null, E: null, extinta: estado.extinta, aviso: null
    })];
    var actual = estado;
    for (var i = 0; i < n; i++) {
      var r = paso(actual);
      actual = r.estado;
      trayectoria.push(r.registro);
    }
    return { parametros: estado.parametros, estadoFinal: actual, trayectoria: trayectoria };
  }

  // ---------------------------------------------------------------------- CSV
  // Devuelve solo el texto. Líneas "# clave,valor" con parámetros y semilla, luego la tabla.
  // opciones: { separador: ',' (o ';'), decimal: '.' (o ',') }
  function exportarCSV(resultado, opciones) {
    if (!resultado || !resultado.parametros || !Array.isArray(resultado.trayectoria)) {
      throw new Error('exportarCSV espera el objeto devuelto por simular().');
    }
    var op = opciones || {};
    var sep = op.separador || ',';
    var dec = op.decimal || '.';
    if (sep === dec) throw new Error('El separador y el símbolo decimal no pueden ser iguales.');
    var p = resultado.parametros;

    function num(x) {
      if (x === null || x === undefined) return '';
      return String(x).replace('.', dec);
    }
    function celda(x) {
      var s = (x === null || x === undefined) ? '' : String(x);
      return /["\n\r]/.test(s) || s.indexOf(sep) >= 0 ? '"' + s.replace(/"/g, '""') + '"' : s;
    }
    var L = [];
    L.push('# modelo' + sep + celda('N(t+1) = N + B - D + I - E'));
    L.push('# modo' + sep + p.modo);
    L.push('# N0' + sep + num(p.N0));
    L.push('# natalidad' + sep + num(p.natalidad));
    L.push('# mortalidad' + sep + num(p.mortalidad));
    L.push('# inmigracion' + sep + num(p.inmigracion));
    L.push('# emigracion' + sep + num(p.emigracion));
    L.push('# densidad' + sep + p.densidad);
    if (p.densidad) L.push('# K' + sep + num(p.K));
    L.push('# semilla' + sep + (p.modo === 'estocastico' ? p.semilla : celda('no usada (modo determinista)')));
    L.push(['t', 'N', 'B', 'D', 'I', 'E', 'extinta', 'aviso'].join(sep));
    resultado.trayectoria.forEach(function (r) {
      L.push([num(r.t), num(r.N), num(r.B), num(r.D), num(r.I), num(r.E), r.extinta, celda(r.aviso)].join(sep));
    });
    return L.join('\n') + '\n';
  }

  // -------------------------------------------------------------- equilibrio
  // Equilibrio teórico del modelo DETERMINISTA (sin umbral de extinción). null si no existe.
  function equilibrio(parametros) {
    var v = validarParametros(parametros);
    if (!v.ok) throw new Error('Parámetros inválidos:\n - ' + v.errores.join('\n - '));
    var p = v.parametros, b = p.natalidad, d = p.mortalidad, e = p.emigracion, I = p.inmigracion;
    var r = b - d - e;
    if (p.densidad && b > 0) {
      var a = b / p.K;                                   // a·N² − r·N − I = 0
      var N = (r + Math.sqrt(r * r + 4 * a * I)) / (2 * a);
      if (!(N > 0)) return null;
      var deriv = 1 + r - 2 * a * N;
      return { N: N, estable: Math.abs(deriv) < 1 };
    }
    if (r < 0) return { N: I / (-r), estable: true };    // sin densidad (o b = 0)
    return null;
  }

  // ------------------------------------------------------------------ pruebas
  function pruebas(opciones) {
    var silencioso = !!(opciones && opciones.silencioso);
    var resultados = [];

    function prueba(nombre, fn) {
      try { resultados.push({ nombre: nombre, ok: true, detalle: fn() || '' }); }
      catch (err) { resultados.push({ nombre: nombre, ok: false, detalle: err.message }); }
    }
    function afirmar(cond, msg) { if (!cond) throw new Error(msg || 'afirmación falsa'); }
    function cercano(a, b, tol, msg) {
      if (!(Math.abs(a - b) <= tol)) throw new Error((msg || '') + ' esperado ≈ ' + b + ', obtenido ' + a + ' (tol ' + tol + ')');
    }
    function correr(p, n) { return simular(crearEstado(p), n).trayectoria; }
    function copiar(base, extra) { return Object.assign({}, base, extra); }
    function media(xs) { return xs.reduce(function (s, x) { return s + x; }, 0) / xs.length; }
    function varianza(xs) { var m = media(xs); return xs.reduce(function (s, x) { return s + (x - m) * (x - m); }, 0) / (xs.length - 1); }
    function lanza(fn) { try { fn(); return null; } catch (err) { return err.message; } }

    prueba('Determinista, b = d, I = E = 0: N constante', function () {
      var tr = correr({ N0: 100, natalidad: 0.07, mortalidad: 0.07 }, 100);
      tr.forEach(function (r) { cercano(r.N, 100, 1e-9, 't=' + r.t); });
      var tr2 = correr({ N0: 100, natalidad: 0.12, mortalidad: 0.05, emigracion: 0.07 }, 100);
      tr2.forEach(function (r) { cercano(r.N, 100, 1e-9, 'con emigración, t=' + r.t); });
      return 'N(100) = ' + tr[100].N;
    });

    prueba('A mano: N=100, b=0,1, d=0,05, I=E=0 → N(1)=105', function () {
      var tr = correr({ N0: 100, natalidad: 0.1, mortalidad: 0.05 }, 1);
      cercano(tr[1].N, 105, 1e-9, 'N(1)');
      cercano(tr[1].B, 10, 1e-9, 'B'); cercano(tr[1].D, 5, 1e-9, 'D');
      afirmar(tr[1].I === 0 && tr[1].E === 0, 'I y E deberían ser 0');
      return 'N(1) = ' + tr[1].N + ', B = ' + tr[1].B + ', D = ' + tr[1].D;
    });

    prueba('Sin natalidad ni inmigración, N nunca crece (ambos modos)', function () {
      var casos = [{ d: 0.1, e: 0 }, { d: 0.2, e: 0.3 }, { d: 0, e: 0.5 }, { d: 0, e: 0 }, { d: 1, e: 0 }];
      var cuenta = 0;
      ['determinista', 'estocastico'].forEach(function (modo) {
        casos.forEach(function (c) {
          for (var s = 1; s <= 10; s++) {
            var tr = correr({ modo: modo, N0: 500, natalidad: 0, mortalidad: c.d, emigracion: c.e, semilla: s }, 40);
            for (var i = 1; i < tr.length; i++) { afirmar(tr[i].N <= tr[i - 1].N, modo + ' d=' + c.d + ' e=' + c.e + ': N creció en t=' + i); cuenta++; }
          }
        });
      });
      return cuenta + ' pasos verificados';
    });

    prueba('Misma semilla → misma trayectoria; semillas distintas → distinta', function () {
      var base = { modo: 'estocastico', N0: 200, natalidad: 0.2, mortalidad: 0.1, inmigracion: 5, emigracion: 0.03 };
      var a = JSON.stringify(correr(copiar(base, { semilla: 42 }), 60));
      var b = JSON.stringify(correr(copiar(base, { semilla: 42 }), 60));
      var c = JSON.stringify(correr(copiar(base, { semilla: 43 }), 60));
      afirmar(a === b, 'misma semilla dio trayectorias distintas');
      afirmar(a !== c, 'semillas distintas dieron la misma trayectoria');
      var e0 = crearEstado(copiar(base, { semilla: 7 }));
      afirmar(JSON.stringify(paso(e0)) === JSON.stringify(paso(e0)), 'paso() no es puro');
      var csvA = exportarCSV(simular(crearEstado(copiar(base, { semilla: 42 })), 30));
      var csvB = exportarCSV(simular(crearEstado(copiar(base, { semilla: 42 })), 30));
      afirmar(csvA === csvB, 'CSV distinto con misma semilla');
      return 'semillas 42 y 43 difieren; 42 y 42 coinciden';
    });

    prueba('N nunca negativa, flujos coherentes (mortalidad/emigración extremas)', function () {
      var casos = [{ d: 1, e: 0 }, { d: 0, e: 1 }, { d: 0.5, e: 0.5 }, { d: 0.99, e: 0.01 }, { d: 0.7, e: 0.3 }];
      var N0s = [1, 7, 1000], cuenta = 0;
      ['determinista', 'estocastico'].forEach(function (modo) {
        casos.forEach(function (c) {
          N0s.forEach(function (N0) {
            for (var s = 1; s <= 20; s++) {
              var tr = correr({ modo: modo, N0: N0, natalidad: 0.3, mortalidad: c.d, emigracion: c.e, inmigracion: 0, semilla: s }, 30);
              for (var i = 1; i < tr.length; i++) {
                var r = tr[i], previo = tr[i - 1].N;
                afirmar(r.N >= 0, 'N negativa');
                afirmar(r.B >= 0 && r.D >= 0 && r.I >= 0 && r.E >= 0, 'flujo negativo');
                afirmar(r.D + r.E <= previo * (1 + 1e-12), 'D + E superó a N(t) (' + modo + ')');
                if (modo === 'estocastico') afirmar(Number.isInteger(r.N), 'N no entera');
                cuenta++;
              }
            }
          });
        });
      });
      return cuenta + ' pasos verificados';
    });

    prueba('Balance de cada paso: N(t+1) = N(t) + B − D + I − E', function () {
      var cuenta = 0;
      ['determinista', 'estocastico'].forEach(function (modo) {
        var tr = correr({ modo: modo, N0: 300, natalidad: 0.25, mortalidad: 0.1, emigracion: 0.05, inmigracion: 7, densidad: true, K: 800, semilla: 9 }, 100);
        for (var i = 1; i < tr.length; i++) {
          var r = tr[i];
          cercano(r.N, tr[i - 1].N + r.B - r.D + r.I - r.E, 1e-9 * Math.max(1, r.N), modo + ' t=' + i);
          cuenta++;
        }
      });
      return cuenta + ' pasos verificados';
    });

    prueba('Densidad (determinista): N → K·(1 − (d+e)/b) sin inmigración', function () {
      var p = { N0: 10, natalidad: 0.3, mortalidad: 0.1, emigracion: 0.05, densidad: true, K: 1000 };
      var esperado = 1000 * (1 - 0.15 / 0.3);                 // 500, calculado a mano
      var tr = correr(p, 500);
      cercano(tr[500].N, esperado, 1e-6, 'desde N0=10');
      cercano(correr(copiar(p, { N0: 3000 }), 500)[500].N, esperado, 1e-6, 'desde N0=3000 > K');
      cercano(equilibrio(p).N, esperado, 1e-9, 'equilibrio()');
      for (var i = 0; i < tr.length; i++) afirmar(tr[i].N >= 0 && tr[i].B >= 0, 'valor negativo');
      return 'N*=' + esperado + ', obtenido ' + tr[500].N.toFixed(8);
    });

    prueba('Densidad (determinista) con I y E: N → raíz positiva de la cuadrática', function () {
      var p = { N0: 50, natalidad: 0.4, mortalidad: 0.1, emigracion: 0.05, inmigracion: 20, densidad: true, K: 1000 };
      var r = 0.25, a = 0.4 / 1000;
      var esperado = (r + Math.sqrt(r * r + 4 * a * 20)) / (2 * a);   // ≈ 696,3
      var tr = correr(p, 600);
      cercano(tr[600].N, esperado, 1e-6, 'N final');
      cercano(equilibrio(p).N, esperado, 1e-9, 'equilibrio()');
      afirmar(equilibrio(p).estable === true, 'debería ser estable');
      return 'N*≈' + esperado.toFixed(4) + ', obtenido ' + tr[600].N.toFixed(4);
    });

    prueba('Densidad (estocástico): N fluctúa alrededor del equilibrio (±3 %)', function () {
      var p = { modo: 'estocastico', N0: 50, natalidad: 0.3, mortalidad: 0.1, emigracion: 0.05, densidad: true, K: 1000 };
      var medias = [];
      for (var s = 1; s <= 20; s++) {
        var tr = correr(copiar(p, { semilla: s }), 600).slice(501);
        medias.push(media(tr.map(function (r) { return r.N; })));
      }
      var m = media(medias);
      afirmar(Math.abs(m - 500) <= 15, 'media ' + m.toFixed(1) + ' lejos de 500');
      return 'media de 20 réplicas (t=500..600) = ' + m.toFixed(1) + ' (equilibrio 500)';
    });

    prueba('Estocástico: la media de muchas réplicas ≈ trayectoria determinista', function () {
      var base = { N0: 1000, natalidad: 0.12, mortalidad: 0.05, emigracion: 0.02, inmigracion: 10 };
      var n = 25, R = 500;
      var det = correr(copiar(base, { modo: 'determinista' }), n)[n].N;
      var finales = [];
      for (var s = 1; s <= R; s++) finales.push(correr(copiar(base, { modo: 'estocastico', semilla: s }), n)[n].N);
      var m = media(finales), ee = Math.sqrt(varianza(finales) / R);
      afirmar(ee > 0, 'sin variabilidad: no es estocástico');
      afirmar(Math.abs(m - det) <= 4 * ee, 'media ' + m.toFixed(2) + ' vs determinista ' + det.toFixed(2) + ' (EE ' + ee.toFixed(2) + ')');
      return 'determinista = ' + det.toFixed(2) + '; media de ' + R + ' réplicas = ' + m.toFixed(2) + ' (EE ' + ee.toFixed(2) + ')';
    });

    prueba('Extinción terminal: no revive, ni por inmigración', function () {
      var tr = correr({ N0: 3, natalidad: 0, mortalidad: 1, inmigracion: 0 }, 5);
      afirmar(tr[1].N === 0 && tr[1].extinta === true, 'debió extinguirse en t=1');
      tr.slice(1).forEach(function (r) { afirmar(r.N === 0, 'N no es 0 en t=' + r.t); });
      ['determinista', 'estocastico'].forEach(function (modo) {
        var t2 = correr({ modo: modo, N0: 0, natalidad: 0.5, mortalidad: 0.1, inmigracion: 10 }, 10);
        t2.forEach(function (r) { afirmar(r.N === 0 && r.extinta === true, modo + ': revivió en t=' + r.t); });
        t2.slice(1).forEach(function (r) { afirmar(r.I === 0 && r.B === 0, 'flujos no nulos tras extinción'); });
      });
      var t3 = correr({ N0: 1.5, natalidad: 0, mortalidad: 0.5 }, 3);   // 1,5 → 0,75 < 1 ⇒ extinta
      afirmar(t3[1].N === 0 && t3[1].aviso !== null, 'umbral determinista');
      return 'N0=0 con I=10 permanece en 0';
    });

    prueba('Validación: parámetros inválidos rechazados con mensaje', function () {
      var base = { N0: 100, natalidad: 0.1, mortalidad: 0.05 };
      var malos = [
        [copiar(base, { N0: -1 }), 'N0'], [copiar(base, { N0: NaN }), 'N0'], [copiar(base, { N0: '100' }), 'N0'],
        [copiar(base, { N0: 0.5 }), 'N0'], [copiar(base, { N0: 2 * Math.pow(2, 53) }), 'N0'],
        [copiar(base, { natalidad: -0.1 }), 'natalidad'], [copiar(base, { natalidad: Infinity }), 'natalidad'],
        [copiar(base, { mortalidad: 1.5 }), 'mortalidad'], [copiar(base, { mortalidad: -0.2 }), 'mortalidad'],
        [copiar(base, { emigracion: 2 }), 'emigracion'], [copiar(base, { inmigracion: -3 }), 'inmigracion'],
        [copiar(base, { mortalidad: 0.7, emigracion: 0.31 }), 'no puede superar 1'],
        [copiar(base, { densidad: true }), 'K'], [copiar(base, { densidad: true, K: 0 }), 'K'],
        [copiar(base, { densidad: 'si' }), 'densidad'], [copiar(base, { modo: 'azar' }), 'modo'],
        [copiar(base, { semilla: NaN }), 'semilla'], [{ natalidad: 0.1, mortalidad: 0.1 }, 'N0'],
        [{ N0: 5, mortalidad: 0.1 }, 'natalidad'], [null, 'objeto']
      ];
      malos.forEach(function (c) {
        var v = validarParametros(c[0]);
        afirmar(v.ok === false, 'aceptó: ' + JSON.stringify(c[0]));
        afirmar(v.errores.join(' ').indexOf(c[1]) >= 0, 'mensaje sin "' + c[1] + '": ' + v.errores.join(' | '));
        afirmar(lanza(function () { crearEstado(c[0]); }) !== null, 'crearEstado no lanzó');
      });
      afirmar(validarParametros(copiar(base, { mortalidad: 0.7, emigracion: 0.3 })).ok, 'd+e=1 debe aceptarse');
      var corr = validarParametros(copiar(base, { modo: 'estocastico', N0: 10.4, semilla: 3.7, K: 50, azar: 1 }));
      afirmar(corr.ok && corr.parametros.N0 === 10 && corr.parametros.semilla === 3, 'corrección de N0/semilla');
      afirmar(corr.advertencias.length === 4, 'se esperaban 4 advertencias: ' + corr.advertencias.join(' | '));
      return malos.length + ' casos inválidos rechazados; 4 correcciones con advertencia';
    });

    prueba('Pureza: no muta el estado ni depende del objeto de parámetros original', function () {
      var p = { modo: 'estocastico', N0: 100, natalidad: 0.1, mortalidad: 0.05, semilla: 5 };
      var e0 = crearEstado(p);
      var antes = JSON.stringify(e0);
      p.natalidad = 9; p.N0 = 1;
      paso(e0); simular(e0, 10);
      afirmar(JSON.stringify(e0) === antes, 'el estado cambió');
      afirmar(e0.parametros.natalidad === 0.1, 'el estado depende del objeto original');
      afirmar(Object.isFrozen(e0), 'estado no congelado');
      return 'estado inmutable';
    });

    prueba('N muy grande: sin pérdida de precisión ni bloqueo; saturación avisada', function () {
      var t0 = Date.now();
      var p = { modo: 'estocastico', N0: 1e12, natalidad: 0.05, mortalidad: 0.04, inmigracion: 1e6, semilla: 11 };
      var tr = correr(p, 100);
      tr.forEach(function (r) { afirmar(Number.isInteger(r.N) && r.N >= 0 && r.N <= LIMITE_N, 'N inválida en t=' + r.t); });
      var det = correr(copiar(p, { modo: 'determinista' }), 100)[100].N;
      var rel = Math.abs(tr[100].N - det) / det;
      afirmar(rel < 1e-4, 'desvío relativo ' + rel);
      var sat = correr({ modo: 'estocastico', N0: 9e15, natalidad: 1, mortalidad: 0, semilla: 2 }, 3);
      afirmar(sat[1].N === LIMITE_N && sat[1].aviso !== null, 'sin saturación avisada');
      var satD = correr({ N0: 10, natalidad: 50, mortalidad: 0 }, 400);
      afirmar(satD.every(function (r) { return isFinite(r.N); }), 'determinista llegó a Infinity');
      var ms = Date.now() - t0;
      afirmar(ms < 5000, 'demasiado lento: ' + ms + ' ms');
      return 'N(100) = ' + tr[100].N + ' (determinista ' + det.toPrecision(12) + '), ' + ms + ' ms';
    });

    prueba('Muestreadores: media y varianza de binomial y Poisson', function () {
      var rng = crearRng(2024), R = 20000, info = [];
      function revisar(nombre, gen, mediaT, varT) {
        var xs = []; for (var i = 0; i < R; i++) xs.push(gen());
        var m = media(xs), v = varianza(xs);
        afirmar(Math.abs(m - mediaT) <= 4 * Math.sqrt(varT / R), nombre + ': media ' + m + ' vs ' + mediaT);
        afirmar(Math.abs(v / varT - 1) <= 0.1, nombre + ': varianza ' + v + ' vs ' + varT);
        info.push(nombre);
      }
      revisar('Bin(20,0.3)', function () { return binomial(rng, 20, 0.3); }, 6, 4.2);
      revisar('Bin(1e6,1e-5)', function () { return binomial(rng, 1e6, 1e-5); }, 10, 10 * (1 - 1e-5));
      revisar('Bin(100,0.9)', function () { return binomial(rng, 100, 0.9); }, 90, 9);
      revisar('Bin(1e9,0.4)', function () { return binomial(rng, 1e9, 0.4); }, 4e8, 1e9 * 0.4 * 0.6);
      revisar('Poisson(3)', function () { return poisson(rng, 3); }, 3, 3);
      revisar('Poisson(25)', function () { return poisson(rng, 25); }, 25, 25);
      revisar('Poisson(1e6)', function () { return poisson(rng, 1e6); }, 1e6, 1e6);
      afirmar(binomial(rng, 0, 0.5) === 0 && binomial(rng, 9, 0) === 0 && binomial(rng, 9, 1) === 9, 'bordes binomial');
      afirmar(poisson(rng, 0) === 0, 'borde Poisson');
      return info.join(', ');
    });

    prueba('Exportación CSV: estructura, semilla y parámetros', function () {
      var p = { modo: 'estocastico', N0: 50, natalidad: 0.2, mortalidad: 0.1, densidad: true, K: 300, semilla: 42 };
      var res = simular(crearEstado(p), 20);
      var csv = exportarCSV(res);
      var lineas = csv.trim().split('\n');
      var iCab = lineas.findIndex(function (l) { return l.indexOf('t,N,B,D,I,E') === 0; });
      afirmar(iCab > 0, 'sin cabecera');
      afirmar(lineas.length - iCab - 1 === 21, 'filas de datos: ' + (lineas.length - iCab - 1));
      afirmar(csv.indexOf('# semilla,42') >= 0 && csv.indexOf('# K,300') >= 0 && csv.indexOf('# modo,estocastico') >= 0, 'faltan metadatos');
      var fila5 = lineas[iCab + 6].split(',');
      afirmar(Number(fila5[0]) === 5 && Number(fila5[1]) === res.trayectoria[5].N, 'fila 5 no coincide');
      var csvD = exportarCSV(simular(crearEstado({ N0: 100, natalidad: 0.1, mortalidad: 0.05 }), 3), { separador: ';', decimal: ',' });
      afirmar(csvD.indexOf('# semilla;no usada') >= 0 && csvD.indexOf('105') >= 0, 'variante con ; y coma decimal');
      afirmar(lanza(function () { exportarCSV({}); }) !== null, 'debe rechazar entrada inválida');
      return lineas.length + ' líneas';
    });

    var pasaron = resultados.filter(function (r) { return r.ok; }).length;
    var salida = { total: resultados.length, pasaron: pasaron, fallaron: resultados.length - pasaron, resultados: resultados };
    if (!silencioso && typeof console !== 'undefined') {
      resultados.forEach(function (r) {
        console.log((r.ok ? '  ✔ ' : '  ✘ ') + r.nombre + (r.detalle ? '\n      → ' + r.detalle : ''));
      });
      console.log('\n' + pasaron + '/' + resultados.length + ' pruebas superadas' + (salida.fallaron ? ' — HAY FALLOS' : ''));
    }
    return salida;
  }

  // --------------------------------------------------------------- exportación
  var API = Object.freeze({
    validarParametros: validarParametros,
    crearEstado: crearEstado,
    paso: paso,
    simular: simular,
    exportarCSV: exportarCSV,
    equilibrio: equilibrio,
    pruebas: pruebas
  });

  if (typeof module !== 'undefined' && module.exports) {
    module.exports = API;
    if (typeof require !== 'undefined' && require.main === module &&
        typeof process !== 'undefined' && process.argv.indexOf('--test') >= 0) {
      var r = pruebas();
      process.exitCode = r.fallaron ? 1 : 0;
    }
  } else {
    raiz.ModeloPoblacional = API;
  }
})(typeof globalThis !== 'undefined' ? globalThis : this);
