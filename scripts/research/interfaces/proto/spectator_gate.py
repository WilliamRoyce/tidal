"""R-1 prototype: a Cobaya Theory that refuses to sample without its derived spectrum.

Never derives, never starts a Wolfram kernel. Couplings are registered MFLike-style
(``mflike/foreground.py:215-225``), priors stay in top-level ``params:``, refusal uses
Cobaya's own ``ComponentNotInstalledError`` idiom (``camb.py:290-296``), and the fingerprint
is reported as the component version so a resume against another spectrum is refused.
"""

import theory_store as ts
from cobaya.component import ComponentNotInstalledError
from cobaya.log import LoggedError
from cobaya.theory import Theory


class SpectatorGate(Theory):
    model: dict
    spectra_path: str | None = None

    @classmethod
    def get_modified_defaults(
        cls, defaults: dict, input_options: dict | None = None
    ) -> dict:
        model = (input_options or {}).get("model") or defaults.get("model") or {}
        defaults["params"] = {
            **(defaults.get("params") or {}),
            **dict.fromkeys(model.get("lagrangian", {})),
        }
        return defaults

    def initialize(self):
        self._theory = ts.theory_from_mapping(self.model)
        self._fp = ts.fingerprint(self._theory)
        entry, stores = ts.find_entry(self._fp, self.spectra_path)
        if entry is None:
            raise ComponentNotInstalledError(
                self.log,
                (
                    f"No derived spectrum for this theory (fingerprint {self._fp[:12]}). Generate it with "
                    "`tidalcosmo derive <your input>.yaml` on a machine with Wolfram (minutes to hours), "
                    f"or copy it from one. Searched: {', '.join(map(str, stores))}"
                ),
            )
        try:
            self._manifest = ts.verify_entry(entry, self._theory)
        except ValueError as err:
            raise LoggedError(
                self.log,
                f"Derived spectrum at {entry} is unusable: {err}. "
                "Re-derive with `tidalcosmo derive <your input>.yaml --force`.",
            ) from err

    def get_can_support_params(self):
        return ts.couplings(self._theory)

    def initialize_with_params(self):
        # Missing couplings are left to Cobaya's own "Requirement ... not provided" error (it fires
        # later); raising here would preempt it. Only an assignment Cobaya does not police is checked.
        unexpected = sorted(set(self.input_params) - set(ts.couplings(self._theory)))
        if unexpected:
            raise LoggedError(
                self.log,
                f"Parameters {unexpected} were assigned but are not couplings of this theory",
            )

    def get_version(self):
        return f"spectrum-{self._fp}"

    def calculate(  # noqa: PLR6301 - Cobaya's hook must be a method
        self,
        state: dict,
        want_derived: bool = True,  # noqa: ARG002 - part of Cobaya's hook signature
        **params: float,
    ) -> None:
        state["couplings"] = dict(params)
