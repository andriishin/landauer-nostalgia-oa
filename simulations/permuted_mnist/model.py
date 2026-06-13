"""ML surrogate «Permuted MNIST» для § 8.4 main.ru.md (статья #2, Informational
Nostalgia). Самодостаточная numpy-реализация.

Что моделируется
----------------
Непрерывное онлайн-обучение MLP (64 -> H -> 10) на потоке примеров рукописных
цифр (sklearn ``load_digits``, мини-аналог MNIST). Входные координаты (пиксели)
переставляются глобальной перестановкой ``pi_t``; на каждом шаге с вероятностью
``rho_p`` перестановка заменяется новой случайной. Метки неизменны. После смены
``pi`` старые веса входного слоя устаревают, и сеть вынуждена переучиваться.

Это инстанцирует канонический параметр статьи
``c := M_refresh_dot * tau_E / |M_{<=t}|`` (§ 6.2/§ 6.3). Для Permuted-задачи:
``tau_E = 1/rho_p`` (характерное время постоянства среды), темп обновления памяти
SGD ``M_refresh_dot ~ eta*N`` бит/шаг, ёмкость ``|M| ~ N`` (N — число параметров)
=> ``c ~ eta/rho_p`` (§ 8.4).

Истинно-ностальгическая фракция nu^op_true
------------------------------------------
ИСПРАВЛЕНО (см. README «Честная рамка»). Прежняя метрика ``nu^op`` = доля юнитов,
чьё зануление не повышает loss на текущей задаче. Эта метрика МЕХАНИЧЕСКИ ->1 при
быстром дрейфе просто потому, что сеть НЕ ВЫУЧИВАЕТ текущую задачу (терять нечего):
она измеряла (1 - выученность), а не ностальгию. Поэтому добавлены два контроля.

1. **Accuracy-gate / dwell.** Ностальгия — это память, которая УДЕРЖИВАЕТСЯ и
   КОГДА-ТО была предсказательной. Это осмысленно ТОЛЬКО когда текущая задача
   реально выучена. Поэтому время пребывания ``dwell ~ 1/rho_p`` должно быть
   >> числа шагов до сходимости на одной перестановке; точки с низкой accuracy
   помечаются ``undertrained`` и исключаются из вывода о ну.

2. **Контроль прошлого (определение истинной ностальгии).** Юнит ностальгичен,
   если его зануление (а) НЕ повышает loss на ТЕКУЩЕЙ перестановке выше порога
   (бесполезен сейчас) И (б) ПОВЫШАЕТ loss на хотя бы одной НЕДАВНЕЙ ПРОШЛОЙ
   перестановке (был полезен раньше — кодирует устаревшую структуру). Это
   отделяет stale-retained от never-useful/dead/redundant.

      nu^op_true = mean[ (dL_current <= thr) & (max_k dL_past_k > thr) ]
      dead_frac  = mean[ (dL_current <= thr) & (max_k dL_past_k <= thr) ]
      useful_frac= mean[ dL_current > thr ]
"""

from __future__ import annotations

import numpy as np


# --------------------------------------------------------------------------- #
# Данные                                                                       #
# --------------------------------------------------------------------------- #

def load_dataset(rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray, int, int, str]:
    """Вернуть (X, y, n_features, n_classes, source).

    Предпочтительно реальные рукописные цифры 8x8 из sklearn (без скачивания).
    Фолбэк — синтетический teacher-датасет с фиксированной линейной разметкой,
    если sklearn недоступен.
    """
    try:
        from sklearn.datasets import load_digits  # локальные данные, без сети
        d = load_digits()
        X = d.data.astype(np.float64)
        y = d.target.astype(np.int64)
        # нормировка пикселей в [0,1] (значения 0..16)
        X = X / 16.0
        return X, y, X.shape[1], 10, "sklearn.load_digits (8x8 рукописные цифры)"
    except Exception:
        # Синтетический teacher: argmax(W0 @ x) над N(0,1)-входами.
        n_features, n_classes, n_samples = 64, 10, 2000
        W0 = rng.normal(0.0, 1.0, size=(n_classes, n_features))
        X = rng.normal(0.0, 1.0, size=(n_samples, n_features))
        y = np.argmax(X @ W0.T, axis=1).astype(np.int64)
        # центрируем входы в [0,1]-подобный масштаб
        X = (X - X.min(0)) / (X.max(0) - X.min(0) + 1e-9)
        return X, y, n_features, n_classes, "synthetic teacher argmax(W0@x)"


# --------------------------------------------------------------------------- #
# MLP: 64 -> H -> 10, softmax + cross-entropy                                  #
# --------------------------------------------------------------------------- #

class MLP:
    """Двухслойный перцептрон в чистом numpy с tanh-скрытым слоем.

    Параметры:
      W1: (H, D), b1: (H,)        входной слой
      W2: (C, H), b2: (C,)        выходной слой (исходящие веса юнитов = столбцы)
    """

    def __init__(self, n_in: int, n_hidden: int, n_out: int,
                 rng: np.random.Generator) -> None:
        self.D, self.H, self.C = n_in, n_hidden, n_out
        # He/Xavier-подобная инициализация
        self.W1 = rng.normal(0.0, np.sqrt(1.0 / n_in), size=(n_hidden, n_in))
        self.b1 = np.zeros(n_hidden)
        self.W2 = rng.normal(0.0, np.sqrt(1.0 / n_hidden), size=(n_out, n_hidden))
        self.b2 = np.zeros(n_out)

    # -- forward --------------------------------------------------------- #
    def forward(self, X: np.ndarray) -> tuple[np.ndarray, dict]:
        """X: (B, D). Вернуть (probs (B,C), cache)."""
        z1 = X @ self.W1.T + self.b1        # (B, H)
        a1 = np.tanh(z1)                    # (B, H)
        z2 = a1 @ self.W2.T + self.b2       # (B, C)
        z2 = z2 - z2.max(axis=1, keepdims=True)
        e = np.exp(z2)
        probs = e / e.sum(axis=1, keepdims=True)
        return probs, {"X": X, "a1": a1}

    # -- loss ------------------------------------------------------------ #
    def cross_entropy(self, X: np.ndarray, y: np.ndarray,
                      eps: float = 1e-12) -> float:
        probs, _ = self.forward(X)
        n = X.shape[0]
        return float(-np.log(np.maximum(probs[np.arange(n), y], eps)).mean())

    def accuracy(self, X: np.ndarray, y: np.ndarray) -> float:
        probs, _ = self.forward(X)
        return float((probs.argmax(axis=1) == y).mean())

    # -- one SGD step on a minibatch ------------------------------------- #
    def sgd_step(self, X: np.ndarray, y: np.ndarray, eta: float) -> None:
        B = X.shape[0]
        probs, cache = self.forward(X)
        a1 = cache["a1"]
        # d L / d z2
        dz2 = probs.copy()
        dz2[np.arange(B), y] -= 1.0
        dz2 /= B                            # (B, C)
        dW2 = dz2.T @ a1                     # (C, H)
        db2 = dz2.sum(axis=0)               # (C,)
        da1 = dz2 @ self.W2                  # (B, H)
        dz1 = da1 * (1.0 - a1 ** 2)          # tanh'
        dW1 = dz1.T @ X                      # (H, D)
        db1 = dz1.sum(axis=0)               # (H,)
        self.W2 -= eta * dW2
        self.b2 -= eta * db2
        self.W1 -= eta * dW1
        self.b1 -= eta * db1


# --------------------------------------------------------------------------- #
# counterfactual ablation -> nu^op_true (контроль прошлого)                    #
# --------------------------------------------------------------------------- #

def ablation_deltas(mlp: MLP, X: np.ndarray, y: np.ndarray) -> tuple[float, np.ndarray]:
    """Для каждого скрытого юнита прирост cross-entropy при занулении его
    исходящих весов (столбца W2[:, h]). Вернуть (L_base, dL (H,))."""
    L_base = mlp.cross_entropy(X, y)
    H = mlp.H
    dL = np.zeros(H)
    saved_col = np.empty(mlp.C)
    for h in range(H):
        saved_col[:] = mlp.W2[:, h]
        mlp.W2[:, h] = 0.0
        dL[h] = mlp.cross_entropy(X, y) - L_base
        mlp.W2[:, h] = saved_col
    return L_base, dL


def true_nostalgic_fraction(
    mlp: MLP,
    X_cur: np.ndarray,
    y_cur: np.ndarray,
    past_batches: list[tuple[np.ndarray, np.ndarray]],
    threshold: float,
) -> dict:
    """Истинно-ностальгическая фракция с контролем прошлого.

    Юнит h:
      * useful   : dL_current[h] > threshold              (полезен сейчас);
      * nostalgic: dL_current[h] <= threshold И max_k dL_past_k[h] > threshold
                   (бесполезен сейчас, но был полезен на недавней прошлой
                   перестановке — кодирует устаревшую структуру);
      * dead     : dL_current[h] <= threshold И max_k dL_past_k[h] <= threshold
                   (мёртв/избыточен — нигде не несёт структуры).

    past_batches: список (X_perm_k, y) для нескольких НЕДАВНИХ прошлых
    перестановок; пиксели уже переставлены соответствующими pi. Если список пуст
    (ещё не было ни одной смены), ностальгия не определена -> вернуть NaN-маркеры.

    Вернуть dict(nu_true, dead_frac, useful_frac, dL_current_mean, dL_past_mean).
    """
    _, dL_cur = ablation_deltas(mlp, X_cur, y_cur)

    if not past_batches:
        nan = float("nan")
        return dict(nu_true=nan, dead_frac=nan,
                    useful_frac=float((dL_cur > threshold).mean()),
                    dL_current_mean=float(dL_cur.mean()),
                    dL_past_mean=nan, defined=False)

    H = mlp.H
    dL_past_max = np.full(H, -np.inf)
    for Xp, yp in past_batches:
        _, dLp = ablation_deltas(mlp, Xp, yp)
        dL_past_max = np.maximum(dL_past_max, dLp)

    useful = dL_cur > threshold
    not_useful = ~useful
    nostalgic = not_useful & (dL_past_max > threshold)
    dead = not_useful & (dL_past_max <= threshold)

    return dict(
        nu_true=float(nostalgic.mean()),
        dead_frac=float(dead.mean()),
        useful_frac=float(useful.mean()),
        dL_current_mean=float(dL_cur.mean()),
        dL_past_mean=float(dL_past_max.mean()),
        defined=True,
    )


# --------------------------------------------------------------------------- #
# Дрейфовая онлайн-симуляция для одного rho_p                                  #
# --------------------------------------------------------------------------- #

def simulate_rho(
    X: np.ndarray,
    y: np.ndarray,
    n_classes: int,
    rho_p: float,
    eta: float,
    n_hidden: int,
    n_steps: int,
    batch_size: int,
    abl_threshold: float,
    measure_every: int,
    abl_sample: int,
    n_past: int,
    rng: np.random.Generator,
) -> dict:
    """Один прогон онлайн-SGD под дрейфом перестановок при заданном rho_p.

    На каждом шаге: с вероятностью rho_p берём новую случайную перестановку
    координат pi (предыдущую кладём в очередь недавних прошлых перестановок,
    длиной n_past); формируем минибатч permuted-входов; делаем SGD-шаг.
    Периодически измеряем (а) accuracy на текущей перестановке (accuracy-gate),
    (б) истинно-ностальгическую фракцию nu^op_true с контролем прошлого.

    Усреднение — по последней половине измерений (после установления режима) и
    только по тем измерениям, где ностальгия определена (>=1 прошлая перестановка
    уже произошла). accuracy усредняется по всем измерениям последней половины.
    """
    D = X.shape[1]
    n_samples = X.shape[0]
    mlp = MLP(D, n_hidden, n_classes, rng)

    pi = rng.permutation(D)
    recent: list[np.ndarray] = []   # недавние прошлые перестановки (свежайшая первой)
    n_switches = 0

    times: list[int] = []
    acc_hist: list[float] = []
    nu_hist: list[float] = []
    dead_hist: list[float] = []
    useful_hist: list[float] = []

    for t in range(1, n_steps + 1):
        if rng.random() < rho_p:
            recent.insert(0, pi.copy())
            recent = recent[:n_past]
            pi = rng.permutation(D)
            n_switches += 1
        idx = rng.integers(0, n_samples, size=batch_size)
        Xb = X[idx][:, pi]
        yb = y[idx]
        mlp.sgd_step(Xb, yb, eta)

        if t % measure_every == 0:
            sidx = rng.integers(0, n_samples, size=min(abl_sample, n_samples))
            Xm = X[sidx][:, pi]
            ym = y[sidx]
            past_batches = [(X[sidx][:, pp], y[sidx]) for pp in recent]
            res = true_nostalgic_fraction(mlp, Xm, ym, past_batches, abl_threshold)
            times.append(t)
            acc_hist.append(mlp.accuracy(Xm, ym))
            nu_hist.append(res["nu_true"])
            dead_hist.append(res["dead_frac"])
            useful_hist.append(res["useful_frac"])

    times_a = np.asarray(times)
    acc_a = np.asarray(acc_hist)
    nu_a = np.asarray(nu_hist)
    dead_a = np.asarray(dead_hist)
    useful_a = np.asarray(useful_hist)

    # последняя половина измерений (установившийся режим)
    half = len(times_a) // 2
    sl = slice(half, None)
    acc_late = float(acc_a[sl].mean())

    nu_late_vals = nu_a[sl][np.isfinite(nu_a[sl])]
    dead_late_vals = dead_a[sl][np.isfinite(dead_a[sl])]
    useful_late_vals = useful_a[sl][np.isfinite(useful_a[sl])]

    def _mean(v):
        return float(v.mean()) if v.size else float("nan")

    def _std(v):
        return float(v.std()) if v.size else float("nan")

    return dict(
        rho_p=rho_p,
        eta=eta,
        c_theory=eta / rho_p,
        dwell=1.0 / rho_p,
        accuracy=acc_late,
        nu_true=_mean(nu_late_vals),
        nu_true_std=_std(nu_late_vals),
        dead_frac=_mean(dead_late_vals),
        useful_frac=_mean(useful_late_vals),
        n_nu_measurements=int(nu_late_vals.size),
        n_switches=n_switches,
        times=times_a,
        acc_hist=acc_a,
        nu_hist=nu_a,
        dead_hist=dead_a,
        useful_hist=useful_a,
    )
